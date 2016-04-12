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
        self.topology.ligand = self.u.select_atoms("resname LDP")
        self.topology.ligand_no_H=self.u.select_atoms("resname LDP and not name H*")
        self.topology.define_ligand(self.topology.ligand)
        self.topology.find_res_to_plot(3.5)
        self.hbonds = HBonds(self.topology,PDB,None,self.topology.ligand,0,30,tests=True)
    def tearDown(self):
        del self.topology
        del self.hbonds
    def test_finding_donors_acceptors(self):
    	self.hbonds.find_donors_and_acceptors_in_ligand()
    	acceptors = ['O1', 'O2', 'N1']
    	assert_equal(acceptors,self.hbonds.acceptors)
    def test_hbonds_for_drawing(self):
    	self.hbonds.find_donors_and_acceptors_in_ligand()
    	self.hbonds.analyse_hbonds(PDB, None, self.topology.ligand,0,30,3.5)
    	for_drawing =[(u'OG', u'LDP708'),
 (u'O', u'LDP708'),
 (u'OD1', u'LDP708'),
 (u'OD2', u'LDP708')]
 		assert_equal(for_drawing,self.hbonds.hbonds_for_drawing)

