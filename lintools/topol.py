"""Module that deals with input data from user. Contains class Topol_Data that loads the topology and trajectory and Config class (NYI) which creates 
and/or reads configuration file.
Domicevica,Paramo and Newport, 2015."""

import MDAnalysis
from MDAnalysis.analysis import distances
import time
import numpy as np
import operator
from utils import pdb2mol2
from rdkit import Chem
import sys



class Topol_Data(object):
    def __init__(self, topology, trajectory=None, ligand_name=None, offset=0):
        self.universe = None
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
            self.universe = MDAnalysis.Universe(topology)
            self.frame_count = 1
        else:
            self.universe = MDAnalysis.Universe(topology, trajectory)
            self.frame_count = self.universe.trajectory.n_frames
    def define_ligand(self,ligand_name):
        self.ligand = ligand_name
        self.ligand.resnames = "LIG"
        self.ligand.resname = "LIG"
        self.ligand.write(str("LIG.pdb"))
        self.pdb = "LIG.pdb"
        pdb2mol2.pdb2mol2(self.pdb)
        self.mol2_file = "LIG_test.mol2"
        self.load_mol2_in_rdkit()
    def load_mol2_in_rdkit(self):
        try:
            self.mol2 = Chem.MolFromMol2File(self.mol2_file,removeHs=False)
            mol = Chem.MolFromSmarts('[$([N;!H0;v3]),$([N;!H0;+1;v4]),$([O,S;H1;+0]),$([n;H1;+0])]')
            self.mol2.GetSubstructMatches(mol, uniquify=1)
        except AttributeError:
            self.mol2 = Chem.MolFromMol2File(self.mol2_file,removeHs=False,sanitize=False)
            self.mol2.UpdatePropertyCache(strict=False)
        if self.mol2 == None:
            print "Exiting. No mol2 file was supplied."
            sys.exit()
            # Kind of a debug
    def renumber_system(self, offset=0):
        self.protein = self.universe.select_atoms("protein")
        self.protein.set_resids(self.protein.resids+int(offset))
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
        self.ligand_no_H = self.ligand.select_atoms("not name H*")
        lig_pos = self.ligand_no_H.positions
        for residue in self.dict_of_plotted_res:
            residue_select= self.universe.select_atoms("resid "+str(self.dict_of_plotted_res[residue][0]))
            res_pos = residue_select.positions
            dist_array = MDAnalysis.analysis.distances.distance_array(lig_pos, res_pos)
            min_values_per_atom={}
            i=-1
            for atom in self.ligand_no_H:
                i+=1
                min_values_per_atom[atom.name]=dist_array[i].min()
            sorted_min_values = sorted(min_values_per_atom.items(), key=operator.itemgetter(1))     
            self.closest_atoms[residue]=[(sorted_min_values[0][0],sorted_min_values[0][1])]     
            if self.hbonds!=None:
                check_hbonds = []
                i=-1
                for res in self.hbonds.hbonds_for_drawing:
                    i+=1
                    if residue==self.hbonds.hbonds_for_drawing[i][1]:
                        check_hbonds.append(self.hbonds.hbonds_for_drawing[i][0])
                if len(check_hbonds) > 0:
                    self.closest_atoms[residue]=[]
                    for atom in check_hbonds:
                        item = atom, min_values_per_atom[atom]
                        self.closest_atoms[residue].append(item)



                            









