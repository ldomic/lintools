from numpy.testing import TestCase, assert_equal, assert_almost_equal
import unittest
import os
from lintools.lintools.topol import Topol_Data
from lintools.lintools.testsuite.datafiles import *
import numpy as np

class TestCheckLigand(TestCase):
    def setUp(self):
        self.filename = PDB
        self.topology = Topol_Data(self.filename)
    def tearDown(self):
        del self.topology
        if os.path.isfile("LIG.pdb")==True:
            os.remove("LIG.pdb")
    def test_renumber(self):
        res1 = self.topology.universe.select_atoms("resid 1")
        self.topology.renumber_system(60)
        res61 = self.topology.universe.select_atoms("resid 61")
        assert_equal(res1.resnames[0],res61.resnames[0])
    def test_find_residues(self):
        self.u = self.topology.universe
        self.topology.ligand = self.u.select_atoms("resname LDP")
        self.topology.find_res_to_plot()
        plotted_res = {'ALA117': [117, 1],
                       'ASP121': [121, 1],
'ASP46': [46, 1],
 'PHE325': [325, 1],
 'PHE43': [43, 1],
 'SER421': [421, 1],
 'SER422': [422, 1],
 'TYR124': [124, 1],
 'VAL120': [120, 1]}
        assert_equal(plotted_res,self.topology.dict_of_plotted_res)
    def test_get_closest_atoms(self):
        self.u = self.topology.universe
        self.topology.ligand = self.u.select_atoms("resname LDP")
        self.topology.find_res_to_plot()
        self.topology.get_closest_ligand_atoms()
        closest_res ={'ASP121': ('O1',2.7923072904423361,'O2',2.9820389652431976,'C3',3.8011160428911479),
 'PHE325': ('C6',
  3.3749625311013798,
  'C5',
  3.6330910060120107,
  'C7',
  4.3617278148294476),
 }
        assert_equal(closest_res["ASP121"],self.topology.closest_atoms["ASP121"])
        assert_equal(closest_res["PHE325"],self.topology.closest_atoms["PHE325"])


