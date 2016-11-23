import os
from rdkit import Chem
from collections import namedtuple, defaultdict
import numpy as np
import MDAnalysis
from timeit import default_timer as timer


class SaltBridges(object):
    """This module analyses salt bridges between protein and ligand. A salt bridge is defined
    as a charged interaction between two charge centers - one on ligand and one on protein.
    The charge centre of ligand is an atom with formal charge that is not equal to 0. 
    The charge centre of protein is on either LYS,ARG,HIS (positive charge) or GLU,ASP (negative charge).
    If the distance between charge centers is lower than cutoff, it is recognised as a salt bridge.
    Takes:
        * topology data object * - information about system (lintools.Data object)
        * trajectory * - list of trajectories
        * ligand_descr_object * - information about ligand properties (lintools.analysis.LigDescr object)
        * start_frame_num * - list^^ of frame numbers for start of analysis (Opt)
        * end_frame_num * - list^^ of frame numbers for end of analysis (Opt) 
        * skip * - list^^ of numbers of how many frames should be skipped for this analysis (Opt)
        * analysis_cutoff * - a fraction of simulation time a feature has to be present for to be
            plotted in the final figure

    ^^ The reason these values are lists is because several trajectories can be submitted for 
    analysis and different values could be required for each simulation. Therefore, similarly 
    as with trajectories, start, end and skip variables are submitted as lists with values 
    corresponding for each trajectory.
    
    Example: trajectory = ["1.xtc","2.xtc"] #Both are 1000 frames, but the user only wants to 
             analyse second half the the second trajectory
             start = [0(for the first traj),500(for the second traj)]   
             Other values can be left as default.
    Output:
        * self.saltbridges * - numpy array of salt bridges searchable by keywords
        * self.saltbridges_by_time * - numpy array of total salt bridges per frame
        * self.saltbridges_by_type * - numpy array of total count of each particular salt bridge
        """
    def __init__(self,topology_data_object,trajectory,ligand_descr_object,start_frame_num,end_frame_num,skip,analysis_cutoff):
        self.lig_descr = ligand_descr_object
        self.topology_data = topology_data_object
        
        self.saltbridge_dist = 5.5
        self.trajectory = trajectory
        self.start_frame_num = start_frame_num
        self.end_frame_num = end_frame_num
        self.skip = skip
        self.saltbridges={}
        self.saltbridges_by_type={}
        self.saltbridges_by_time={}

        self.find_salt_bridges()
        self.get_saltbridge_frequency(analysis_cutoff)
    def find_salt_bridges(self):
        """First the atom names of protein charge centers are found as these can be different in different forcefields. 
        Then a Euclidian distance is measured between charged ligand atoms and respective inverted charge centres on
        protein. If the distance s lower than threshold (5.5 angstrom) the interaction is recognised as a salt bridge."""
        #Define charge centers for protein - atom names 
        charge_definitions = {"ARG":"name C* and around 1.4 (name N* and not name N)", #carbon atom surrounded by nitrogens that are not in the peptide bond
                            "LYS":"name N* and not name N", # The nitrogen that is not in the peptide bond
                            "HIS":"name N* and not name N", # Both nitrogens (charge distributed between them)
                            "ASP":"name C* and around 1.4 (name O* and not name O)", # Carbon surrounded by oxygens not in peptide bond
                            "GLU":"name C* and around 1.4 (name O* and not name O)" # Carbon surrounded by oxygens not in peptide bond
        }
        prot_charge_center = {}
        # Find the atom names - these can different in different forcefields
        for res in ["ASP","GLU","HIS","LYS","ARG"]:
            for residue in self.topology_data.universe.residues:
                if residue.resname ==res:
                    atomselection = residue.atoms
                    name_selection = atomselection.select_atoms(charge_definitions[res])
                    prot_charge_center[res]=name_selection.names[0]
                    break

        #Measure distances
        data = namedtuple("saltbridge","frame time ligandatomid ligandatomname distance resname resid segid")
        i=0
        if self.trajectory==[]:
            self.trajectory = [self.topology_data.universe.filename]
            self.start_frame_num=[None]
            self.end_frame_num = [None]
            self.skip =[None]
        for traj in self.trajectory:
            self.timeseries=[]
            self.timesteps=[frame.time for frame in self.topology_data.universe.trajectory[self.start_frame_num[i]:self.end_frame_num[i]:self.skip[i]]]
            start = timer()
            self.topology_data.load_trajectory(traj)
            for atom in self.lig_descr.ligand_atoms:
                if self.lig_descr.ligand_atoms[atom]["Formal charges"]<0:
                    for residue in self.topology_data.dict_of_plotted_res:
                        if residue[0] in ["LYS","ARG","HIS"]:
                            pos_res_atom = self.topology_data.universe.select_atoms("resname "+residue[0]+" and resid "+residue[1]+" and segid "+ residue[2]+" and name "+prot_charge_center[residue[0]])
                            lig_atom = self.topology_data.universe.ligand.select_atoms("name "+self.lig_descr.ligand_atoms[atom]["name"])
                            for frame in self.topology_data.universe.trajectory[self.start_frame_num[i]:self.end_frame_num[i]:self.skip[i]]:
                                if residue[0]==["HIS"]:
                                    #HIS has two atoms in the selection - need to find a middle point (i.e. centroid)
                                    pos_res_pos = pos_res_atom.centroid()
                                    dist = self.euclidean3d(pos_res_pos,lig_atom.positions[0])
                                else:
                                    dist = self.euclidean3d(pos_res_atom.positions[0],lig_atom.positions[0])
                                if dist <= self.saltbridge_dist:
                                    contacts = data(frame=frame.frame, time=frame.time, ligandatomid=lig_atom.atoms.ids, ligandatomname=lig_atom.atoms.names,
                                                    distance=dist, resname=residue[0],resid=residue[1],segid=residue[2])
                                    self.timeseries.append(contacts)
                if self.lig_descr.ligand_atoms[atom]["Formal charges"]>0:
                    for residue in self.topology_data.dict_of_plotted_res:
                        if residue[0] in ["ASP","GLU"]:
                            neg_res_atom = self.topology_data.universe.select_atoms("resname "+residue[0]+" and resid "+residue[1]+" and segid "+ residue[2]+" and name "+prot_charge_center[residue[0]])
                            lig_atom = self.topology_data.universe.ligand.select_atoms("name "+self.lig_descr.ligand_atoms[atom]["name"])
                            for frame in self.topology_data.universe.trajectory[self.start_frame_num[i]:self.end_frame_num[i]:self.skip[i]]:
                                dist = self.euclidean3d(neg_res_atom.positions[0],lig_atom.positions[0])
                                if dist <= self.saltbridge_dist:
                                    contacts = data(frame=frame.frame, time=frame.time, ligandatomid=lig_atom.atoms.ids[0], ligandatomname=lig_atom.atoms.names[0],
                                                    distance=dist, resname=residue[0],resid=residue[1],segid=residue[2])
                                    self.timeseries.append(contacts)
            self.saltbridges[i] = self.make_table()

            self.saltbridges_by_time[i] = self.count_by_time()
            self.saltbridges_by_type[i] = self.count_by_type()
            i+=1
        end = timer()
        print "Salt Bridges:"+str(end-start)

    def make_table(self):
        """Make numpy array from timeseries data."""
        num_records = np.sum([1 for frame in self.timeseries])
        dtype = [("frame",float),("time",float),("ligand atom id",int),
                ("ligand atom name","|U4"),("distance",float),
                ("resid",int),("resname","|U4"),("segid","|U4") ]
        out = np.empty((num_records,),dtype=dtype)
        cursor=0
        for contact in self.timeseries:
            out[cursor] = (contact.frame, contact.time,contact.ligandatomid,contact.ligandatomname,contact.distance,
                           contact.resid,contact.resname,contact.segid)
            cursor+=1
        return out.view(np.recarray)

    def count_by_type(self):
        """Count how many times each individual salt bridge occured throughout the simulation.
        Returns numpy array."""
        saltbridges = defaultdict(int)
        for contact in self.timeseries:
            #count by residue name not by proteinring
            pkey = (contact.ligandatomid,contact.ligandatomname, contact.resid,contact.resname,contact.segid)
            saltbridges[pkey]+=1
        dtype = [("ligand_atom_id",int),("ligand_atom_name","|U4"),("resid",int),("resname","|U4"),("segid","|U4"),("frequency",float) ]
        out = np.empty((len(saltbridges),),dtype=dtype)
        tsteps = float(len(self.timesteps))
        for cursor,(key,count) in enumerate(saltbridges.iteritems()):
            out[cursor] = key + (count / tsteps,)
        return out.view(np.recarray)

    def count_by_time(self):
        """Count how many salt bridges occured in each frame. 
        Returns numpy array."""
        out = np.empty((len(self.timesteps),), dtype=[('time', float), ('count', int)])
        for cursor,timestep in enumerate(self.timesteps):
            out[cursor] = (timestep,len([x for x in self.timeseries if x.time==timestep]))
        return out.view(np.recarray)

    def get_saltbridge_frequency(self,analysis_cutoff):
        """Calculates the frequency of salt bridges throughout simulations. If the frequency exceeds the 
        analysis cutoff, this interaction will be taken for further consideration.
        Takes:
            * analysis_cutoff * - fraction of simulation time a feature has to be present for to be plotted
        Output:
            * self.saltbridge_frequency * - frequency of each salt bridge """
        self.frequency = defaultdict(int)
        for traj in self.saltbridges_by_type:
            for contact in self.saltbridges_by_type[traj]:
                self.frequency[tuple(map(tuple,contact["ligandatomid"])),contact["ligandatomname"],contact["resid"],contact["resname"],contact["segid"]]+=contact["frequency"]
        draw_frequency = {i:self.frequency[i] for i in self.frequency if self.frequency[i]>(int(len(self.trajectory))*analysis_cutoff)}

        self.saltbridges_for_drawing = {}
        for contact in draw_frequency:
            self.saltbridges_for_drawing[contact]=draw_frequency[contact]

    def euclidean3d(self,v1, v2):
        """Faster implementation of euclidean distance for the 3D case."""
        if not len(v1) == 3 and len(v2) == 3:
            print("Vectors are not in 3D space. Returning None.")
            return None
        return np.sqrt((v1[0] - v2[0]) ** 2 + (v1[1] - v2[1]) ** 2 + (v1[2] - v2[2]) ** 2)