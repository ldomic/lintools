from numpy.testing import TestCase, assert_equal, assert_almost_equal
import unittest
import os
from lintools.lintools.topol import Topol_Data
from lintools.lintools.analysis.occurrence import Occurrence_analysis
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
        self.occurrence = Occurrence_analysis(TOPOLOGY, [TRAJ_20FR,TRAJ_50FR], self.ligand, 3.5, 0, self.topology)
    def tearDown(self):
        del self.topology
        del self.occurrence
    def test_occurrence_analysis(self):
        self.occurrence.get_closest_residues(30)
        assert_equal(len(self.occurrence.residue_counts),2)
        assert_equal(self.occurrence.universe.frame_count,[21, 51])
        assert_equal(self.occurrence.residue_counts[1]["PHE953"],17)
        assert_equal(self.occurrence.residue_counts[2]["PHE953"],43)