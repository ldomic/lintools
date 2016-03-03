"""Module that deals with input data from user. Contains class Topol_Data that loads the topology and trajectory and Config class (NYI) which creates 
and/or reads configuration file.
Domicevica,Paramo and Newport, 2015."""

import MDAnalysis
from MDAnalysis.analysis import distances
import time
import numpy as np
import operator
import openbabel


class Topol_Data(object):
    def __init__(self, topology, trajectory=None, ligand_name=None, offset=0):
        self.universe = None
        self.protein = None
        self.ligand = None
        self.ligand_no_H =None
        self.protein_selection =None
        self.frame_count=None
        self.closest_atoms={}
        self.dict_of_plotted_res={}
        self.load_system(topology, trajectory)
        self.renumber_system(offset)
    def load_system(self, topology, trajectory):
        if trajectory is None:
            self.topology=MDAnalysis.Universe(topology)
            self.universe = MDAnalysis.Universe(topology)
            self.frame_count = self.universe.trajectory.n_frames
        else:
            self.topology=MDAnalysis.Universe(topology)
            self.universe = MDAnalysis.Universe(topology, trajectory)
            self.frame_count = self.universe.trajectory.n_frames
    def define_ligand(self,ligand_name):
        self.ligand = ligand_name
        self.ligand.resnames = "LIG"
        self.ligand.resname = "LIG"
        self.ligand.write(str("LIG.pdb"))
        self.pdb ="LIG.pdb"
    def make_mol2_file(self):
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("pdb","mol2")
        mol = openbabel.OBMol()
        obConversion.ReadFile(mol,"LIG.pdb")
        obConversion.WriteFile(mol, "LIG.mol2")
        self.mol2_file = "LIG.mol2"

    def renumber_system(self, offset=0):
        self.protein = self.universe.select_atoms("protein")
        self.protein.set_resids(self.protein.resids+int(offset))
        protein = self.topology.select_atoms("protein")
        protein.set_resids(protein.resids+int(offset))
    def find_res_to_plot(self, cutoff=3.5):
        self.protein_selection = self.universe.select_atoms('protein and around '+str(cutoff)+' (segid '+str(self.ligand.segids[0])+' and resid '+str(self.ligand.resids[0])+')')
        for atom in self.protein_selection:
            if atom.resid  not in self.dict_of_plotted_res.keys():
                #for non-analysis plots
                self.dict_of_plotted_res[atom.resname+str(atom.resid)]=[atom.resid, 1]
            
    def get_closest_ligand_atoms(self, hbond_object=None):
        """Finds the ligand atom that is closest to a particular residue"""
        ## Selecting ligand without hydrogen atoms as these are not depicted in the RDKit fig
        self.hbonds=hbond_object
        self.ligand_no_H = self.topology.select_atoms('segid '+str(self.ligand.segids[0])+' and resid '+str(self.ligand.resids[0])+" and not name H*")
        lig_pos = self.ligand_no_H.positions
        for residue in self.dict_of_plotted_res:
            residue_select= self.topology.select_atoms("resid "+str(self.dict_of_plotted_res[residue][0]))
            res_pos = residue_select.positions
            dist_array = MDAnalysis.analysis.distances.distance_array(lig_pos, res_pos)
            min_values_per_atom={}
            i=-1
            for atom in self.ligand_no_H:
                i+=1
                min_values_per_atom[atom.name]=dist_array[i].min()

            sorted_min_values = sorted(min_values_per_atom.items(), key=operator.itemgetter(1))     
            self.closest_atoms[residue]=sorted_min_values[0][0],sorted_min_values[0][1],sorted_min_values[1][0],sorted_min_values[1][1],sorted_min_values[2][0],sorted_min_values[2][1]       
        if self.hbonds!=None:
            check_hbonds={}
            for atom in self.closest_atoms:
                check_hbonds[atom]=[]
                i=-1
                for res in self.hbonds.hbonds_for_drawing:
                    i+=1
                    if atom ==self.hbonds.hbonds_for_drawing[i][1]:
                        check_hbonds[atom].append(self.hbonds.hbonds_for_drawing[i][0])
            
            for atom in check_hbonds:
                if len(check_hbonds[atom])==1 and check_hbonds[atom]!=self.closest_atoms[atom][0]:
                    self.closest_atoms[atom]=check_hbonds[atom][0],self.hbonds.distance
                if len(check_hbonds[atom])==2:
                    closest_dist=self.closest_atoms[atom]
                    if check_hbonds[atom][0]==closest_dist[0]:
                        self.closest_atoms[atom]=closest_dist[0],closest_dist[1],check_hbonds[atom][1],self.hbonds.distance
                    else:
                        self.closest_atoms[atom]=closest_dist[0],closest_dist[1],check_hbonds[atom][0],self.hbonds.distance
                if len(check_hbonds[atom])==3:
                    closest_dist=self.closest_atoms[atom]
                    if check_hbonds[atom][0]==closest_dist[0]:
                        self.closest_atoms[atom]=closest_dist[0],closest_dist[1],check_hbonds[atom][1],self.hbonds.distance, check_hbonds[atom][2],self.hbonds.distance
                    if check_hbonds[atom][1]==closest_dist[0]:
                        self.closest_atoms[atom]=closest_dist[0],closest_dist[1],check_hbonds[atom][0],self.hbonds.distance, check_hbonds[atom][2],self.hbonds.distance
                    if check_hbonds[atom][2]==closest_dist[0]:
                        self.closest_atoms[atom]=closest_dist[0],closest_dist[1],check_hbonds[atom][0],self.hbonds.distance, check_hbonds[atom][1],self.hbonds.distance





class Config(object):
    def __init__(self, topol_object=None):
        self.topol_obj=topol_object
        self.topology=None
        self.trajectory=[]
        #self.write_config_file(outname,topology,trajectory, res_offset, molecule_file, diagram_type, cutoff_distance, analysis_type, domain_file)
    def write_config_file(self, outname, topology, trajectory, res_offset, diagram_type, cutoff_distance, analysis_type, domain_file,analysis_cutoff):
        config_file=open(outname+"_config.txt", "w")
        config_file.write("This log file was created at "+time.strftime("%c")+" and saved as "+outname+"_config.txt. \n \n \n")
        config_file.write("Topology input file:  "+topology+"\n \n")
        config_file.write("Trajectory input file(s):  "+str(trajectory)+"\n \n")
        config_file.write("Residue offset:  "+str(res_offset)+" \n\n")
        config_file.write("Selected ligand residue:  "+str(self.topol_obj.ligand.resnames[0])+" "+str(self.topol_obj.ligand.resids[0])+" on chain "+str(self.topol_obj.ligand.segids[0])+"\n \n")
        config_file.write("Diagram type:  "+diagram_type+"\n \n")
        config_file.write("Cutoff distance:  "+str(cutoff_distance)+" Angstrom \n \n")
        config_file.write("Type of analysis used:  "+str(analysis_type)+"\n\n")
        config_file.write("Domain file:  "+str(domain_file)+"\n \n")
        config_file.write("Analysis cutoff:  "+str(analysis_cutoff)+"\n \n")
        config_file.close()
    def read_config_file(self, config_file):
        with open(config_file,"r") as f:
            lines = f.readlines()
            for line in lines:
                if line.startswith("Topology input file:"):
                    self.topology = line.rsplit(":",2)[1][2:-1]
                if line.startswith("Trajectory input file(s):  "):
                    my_line= line.rsplit(":",2)[1][3:-2]
                    for i in range(len(my_line.rsplit(", ",3))):
                        if i==0:
                            self.trajectory.append(my_line.rsplit(",",3)[i][1:-1])
                        else:
                            self.trajectory.append(my_line.rsplit(",",3)[i][2:-1])
                if line.startswith("Residue offset:"):
                    self.res_offset=line.rsplit(":",2)[1][2:-1]
                if line.startswith("Diagram type:"):
                    self.diagram_type=line.rsplit(":",2)[1][2:-1]
                if line.startswith("Cutoff distance:"):
                    self.cutoff=line.rsplit(":",2)[1][2:-11]
                if line.startswith("Type of analysis used:"):
                    self.analysis_type=line.rsplit(":",2)[1][2:-1]
                if line.startswith("Domain file:"):
                    self.domain_file=line.rsplit(":",2)[1][2:-1]
                if line.startswith("Analysis cutoff:"):
                    self.domain_file=line.rsplit(":",2)[1][2:-1]
                if line.startswith("Selected ligand residue"):
                    self.ligand_name=[line.rsplit(":",2)[1][1:].rsplit(" ",5)[1], line.rsplit(":",2)[1][1:].rsplit(" ",5)[2], line.rsplit(":",2)[1][1:-1].rsplit(" ",5)[5]]
                    

                            









