"""Module that deals with input data from user. Contains class Topol_Data that loads the topology and trajectory and Config class (NYI) which creates 
and/or reads configuration file.
Domicevica,Paramo and Newport, 2015."""

import MDAnalysis
from MDAnalysis.analysis import distances
import time



class Topol_Data(object):
    def __init__(self, topology, trajectory=None, ligand_name=None, offset=0):
        self.universe = None
        self.protein = None
        self.ligand = None
        self.ligand_no_H =None
        self.closest_atoms={}
        self.dict_of_plotted_res={}
        self.load_system(topology, trajectory)
        self.define_ligand(ligand_name)
        self.renumber_system(offset)
    def load_system(self, topology, trajectory):
        if trajectory is None:
            self.universe = MDAnalysis.Universe(topology)
        else:
            self.universe = MDAnalysis.Universe(topology, trajectory)
    def define_ligand(self,ligand_name):
        self.ligand = ligand_name
    def renumber_system(self, offset=0):
        self.protein = self.universe.select_atoms("protein")
        self.protein.set_resids(self.protein.resids+int(offset))
    def find_res_to_plot(self, cutoff=3.5):
        selection = self.universe.select_atoms('protein and around '+str(cutoff)+' (segid '+str(self.ligand.segids[0])+' and resid '+str(self.ligand.resids[0])+')')
        for atom in selection:
            if atom.resid not in self.dict_of_plotted_res.keys():
                #for non-analysis plots
                self.dict_of_plotted_res[atom.resname+str(atom.resid)]=[atom.resid, 1]
            
    def get_closest_ligand_atoms(self):
        """Finds the ligand atom that is closest to a particular residue"""
        ## Selecting ligand without hydrogen atoms as these are not depicted in the RDKit fig
        self.ligand_no_H = self.universe.select_atoms('segid '+str(self.ligand.segids[0])+' and resid '+str(self.ligand.resids[0])+" and not name H*")
        lig_pos = self.ligand_no_H.positions
        for residue in self.dict_of_plotted_res:
            residue_select= self.universe.select_atoms("resid "+str(self.dict_of_plotted_res[residue][0]))
            res_pos = residue_select.positions
            dist_array = MDAnalysis.analysis.distances.distance_array(lig_pos, res_pos)
            i=-1
            for atom in self.ligand_no_H:
                i+=1
                if dist_array[atom.id].min()== dist_array.min():
                    self.closest_atoms[residue]=atom.name,i, dist_array[atom.id].min()

class Config(object):
    def __init__(self, topol_object,outname, topology,trajectory, res_offset, molecule_file, diagram_type, cutoff_distance, analysis_type, domain_file):
        self.topol_obj=topol_object
        self.write_config_file(outname,topology,trajectory, res_offset, molecule_file, diagram_type, cutoff_distance, analysis_type, domain_file)
    def write_config_file(self, outname, topology, trajectory, res_offset, molecule_file, diagram_type, cutoff_distance, analysis_type, domain_file):
        config_file=open(outname+"_config.txt", "w")
        #Greetings, command information, etc
        config_file.write("This log file was created at "+time.strftime("%c")+" and saved as "+outname+"_config.txt. \n \n \n")
        config_file.write("Topology input file:  "+topology+"\n \n")
        config_file.write("Trajectory input file(s):  "+str(trajectory)+"\n \n")
        config_file.write("Residue offset: "+str(res_offset)+" \n \n")
        config_file.write("Molecule input file:  "+molecule_file+"\n \n")
        config_file.write("Selected ligand residue : "+str(self.topol_obj.ligand.resnames[0])+" "+str(self.topol_obj.ligand.resids[0])+" on chain "+str(self.topol_obj.ligand.segids[0])+"\n \n")
        config_file.write("Diagram type: "+diagram_type+"\n \n")
        config_file.write("Cutoff distance: "+str(cutoff_distance)+" Angstrom \n \n")
        config_file.write("Type of analysis used: "+str(analysis_type)+" \n \n")
        config_file.write("Domain file: "+str(domain_file)+"\n \n")









