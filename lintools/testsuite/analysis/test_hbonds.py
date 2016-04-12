from numpy.testing import TestCase, assert_equal, assert_almost_equal
import unittest
import os
from lintools.lintools.topol import Topol_Data
from lintools.lintools.analysis.hbonds import HBonds
from lintools.lintools.testsuite.datafiles import *
import numpy as np

class TestHBonds(TestCase):
    def setUp(self):
        self.topology = Topol_Data(PDB,None,None,0)
        self.u = self.topology.universe
        self.ligand = self.u.select_atoms("resname LDP")
        self.ligand.resname = "LIG"
        self.ligand.resnames = "LIG"
        self.topology.define_ligand(self.ligand)
        self.topology.find_res_to_plot(3.5)
        self.hbonds = HBonds(self.topology,PDB,None,self.ligand,0,30,tests=True)
    def tearDown(self):
        del self.topology
        del self.hbonds
    def test_finding_donors_acceptors(self):
    	self.hbonds.find_donors_and_acceptors_in_ligand()
    	acceptors = ['O1', 'O2', 'N1']
    	assert_equal(acceptors,self.hbonds.acceptors)
    def test_hbonds_for_drawing(self):
    	self.hbonds.find_donors_and_acceptors_in_ligand()
    	self.hbonds.analyse_hbonds(PDB, None, self.ligand,0,30,3.5)
    	for_drawing =[('O1', 'SER422'), ('N1', 'ASP46'), ('O2', 'ALA117')]
 	assert_equal(for_drawing,self.hbonds.hbonds_for_drawing)
    def test_get_closest_atoms_with_hbonds(self):
 	    self.hbonds.find_donors_and_acceptors_in_ligand()
    	self.hbonds.analyse_hbonds(PDB, None, self.ligand,0,30,3.5)
    	self.topology.get_closest_ligand_atoms(self.hbonds)
    	closest_atoms = {'ALA117': ('O2', 3.5),
 'ASP121': ('O1', 2.7923072904423361),
 'ASP46': ('N1', 3.5),
 'PHE325': ('C6', 3.3749625311013798),
 'PHE43': ('C7', 4.1442749657414888),
 'SER421': ('C2', 4.0928581319641619),
 'SER422': ('O1', 3.5),
 'TYR124': ('C8', 3.2547764210794483),
 'VAL120': ('C4', 3.5782832437246794)}
    	assert_equal(closest_atoms,self.topology.closest_atoms)


