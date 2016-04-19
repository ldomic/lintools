from numpy.testing import TestCase, assert_equal, assert_almost_equal
import unittest
import os
from lintools.lintools.topol import Topol_Data
from lintools.lintools.analysis.rmsf import RMSF_measurements
from lintools.lintools.testsuite.datafiles import *
import numpy as np

class TestRMSF(TestCase):
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
        ligand_rmsf = {0: 4.1387446317202752,
 1: 4.6187253957045238,
 2: 3.2753842028754034,
 3: 4.1615701267038538,
 4: 3.1328739558388579,
 5: 3.5006571798063768,
 6: 4.7801324647381387,
 7: 4.3159029823823722,
 8: 3.6244073500508804,
 9: 4.4794485431224604,
 10: 3.8519845041292289,
 11: 4.5277155036512893,
 12: 3.7688544007660312,
 13: 3.1940222839172296,
 14: 3.9479333898542839,
 15: 3.4607580891424696,
 16: 3.4922654419157189,
 17: 3.2962539265042383,
 18: 4.3157197648239212,
 19: 4.6448005677759809,
 20: 3.9782854660536255}
        assert_equal(self.rmsf.ligand_rmsf,ligand_rmsf)
        assert_equal(os.path.isfile("rmsf_rmsf_results.dat"),True)