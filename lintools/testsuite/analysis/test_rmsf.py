from numpy.testing import TestCase, assert_equal, assert_almost_equal
import unittest
import os
from lintools.lintools.topol import Topol_Data
from lintools.lintools.analysis.rmsf import RMSF_measurements
from lintools.lintools.testsuite.datafiles import *
import numpy as np

class TestHBonds(TestCase):
    def setUp(self):
        self.topology = Topol_Data(TOPOLOGY,[TRAJ_20FR,TRAJ_50FR],None,0)
        self.u = self.topology.universe
        self.ligand = self.u.select_atoms("resname UNK")
        self.ligand.resname = "LIG"
        self.ligand.resnames = "LIG"
        self.topology.define_ligand(self.ligand)
        self.rmsf = RMSF_measurements(self.topology,TOPOLOGY, [TRAJ_20FR,TRAJ_50FR], self.ligand, 0, "rmsf")
    def tearDown(self):
        del self.topology
        del self.rmsf
    def test_rmsf_measurements(self):
        ligand_rmsf = {0: 0.70174656505967681,
 1: 0.54160793667284157,
 2: 0.68150060667724111,
 3: 0.42721748526707848,
 4: 0.60829760810859757,
 5: 0.60935825062319049,
 6: 0.55657479578842295,
 7: 0.49334593705657837,
 8: 0.453712249216698,
 9: 0.39473035866417011,
 10: 0.31360776926074818,
 11: 0.72792724551009669,
 12: 0.27115441367245063,
 13: 0.43536997263912658,
 14: 0.20199157985533761,
 15: 0.35059617459554393,
 16: 0.48738489351016084,
 17: 0.49290834023913199,
 18: 0.67984972126322585,
 19: 1.5564568419096925,
 20: 1.6772786845332461}
        assert_equal(self.rmsf.ligand_rmsf,ligand_rmsf)
        assert_equal(os.path.isfile("rmsf_rmsf_results.dat"),True)