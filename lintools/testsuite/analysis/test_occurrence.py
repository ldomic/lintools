from numpy.testing import TestCase, assert_equal, assert_almost_equal
import unittest
import os
from lintools.lintools.topol import Topol_Data
from lintools.lintools.analysis.occurrence import Occurrence_analysis
from lintools.lintools.testsuite.datafiles import *
import numpy as np

class TestOccurrence(TestCase):
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
        dict_of_closest_res = {'GLN916': ('916', 7, 17),
 'ILE310': ('310', 17, 38),
 'LEU35': ('35', 21, 49),
 'MET38': ('38', 15, 39),
 'MET39': ('39', 12, 29),
 'MET919': ('919', 17, 44),
 'PHE306': ('306', 18, 48),
 'PHE698': ('698', 11, 24),
 'PHE702': ('702', 8, 17),
 'PHE948': ('948', 12, 29),
 'PHE953': ('953', 17, 43),
 'TYR280': ('280', 13, 28),
 'TYR920': ('920', 7, 18),
 'TYR923': ('923', 15, 42)}
        assert_equal(self.topology.dict_of_plotted_res,dict_of_closest_res)

class TestOccurrence1Traj(TestCase):
    def setUp(self):
        self.topology = Topol_Data(TOPOLOGY,[TRAJ_20FR],None,0)
        self.u = self.topology.universe
        self.ligand = self.u.select_atoms("resname UNK")
        self.ligand.resname = "LIG"
        self.ligand.resnames = "LIG"
        self.topology.define_ligand(self.ligand)
        self.occurrence = Occurrence_analysis(TOPOLOGY, [TRAJ_20FR], self.ligand, 3.5, 0, self.topology)
    def tearDown(self):
        del self.topology
        del self.occurrence
    def test_occurrence_analysis_single_traj(self):
        self.occurrence.get_closest_residues(30)
        assert_equal(len(self.occurrence.residue_counts),1)
        assert_equal(self.occurrence.universe.frame_count,[21])
        assert_equal(self.occurrence.residue_counts[1]["PHE953"],17)
        dict_of_closest_res = {'GLN916': ('916', 7),
 'ILE310': ('310', 17),
 'LEU35': ('35', 21),
 'MET38': ('38', 15),
 'MET39': ('39', 12),
 'MET919': ('919', 17),
 'PHE306': ('306', 18),
 'PHE698': ('698', 11),
 'PHE702': ('702', 8),
 'PHE948': ('948', 12),
 'PHE953': ('953', 17),
 'TYR280': ('280', 13),
 'TYR920': ('920', 7),
 'TYR923': ('923', 15)}
        assert_equal(self.topology.dict_of_plotted_res,dict_of_closest_res)

class TestOccurrence5Traj(TestCase):
    def setUp(self):
        self.topology = Topol_Data(TOPOLOGY,None,None,0)
        self.u = self.topology.universe
        self.ligand = self.u.select_atoms("resname UNK")
        self.ligand.resname = "LIG"
        self.ligand.resnames = "LIG"
        self.topology.define_ligand(self.ligand)
        self.occurrence = Occurrence_analysis(TOPOLOGY, [TRAJ_20FR,TRAJ_50FR,TRAJ_20FR,TRAJ_20FR,TRAJ_50FR], self.ligand, 3.5, 0, self.topology)
    def tearDown(self):
        del self.topology
        del self.occurrence
    def test_occurrence_analysis_single_traj(self):
        self.occurrence.get_closest_residues(30)
        assert_equal(len(self.occurrence.residue_counts),5)
        assert_equal(self.occurrence.universe.frame_count,[21,51,21,21,51])
        assert_equal(self.occurrence.residue_counts[1]["PHE953"],17)
        assert_equal(self.occurrence.residue_counts[2]["PHE953"],43)
        assert_equal(self.occurrence.residue_counts[3]["PHE953"],17)
        assert_equal(self.occurrence.residue_counts[4]["PHE953"],17)
        assert_equal(self.occurrence.residue_counts[5]["PHE953"],43)
        dict_of_closest_res={'GLN916': ('916', 7, 17, 7, 7, 17),
 'ILE310': ('310', 17, 38, 17, 17, 38),
 'LEU35': ('35', 21, 49, 21, 21, 49),
 'MET38': ('38', 15, 39, 15, 15, 39),
 'MET39': ('39', 12, 29, 12, 12, 29),
 'MET919': ('919', 17, 44, 17, 17, 44),
 'PHE306': ('306', 18, 48, 18, 18, 48),
 'PHE698': ('698', 11, 24, 11, 11, 24),
 'PHE702': ('702', 8, 17, 8, 8, 17),
 'PHE948': ('948', 12, 29, 12, 12, 29),
 'PHE953': ('953', 17, 43, 17, 17, 43),
 'TYR280': ('280', 13, 28, 13, 13, 28),
 'TYR920': ('920', 7, 18, 7, 7, 18),
 'TYR923': ('923', 15, 42, 15, 15, 42)}
        assert_equal(self.topology.dict_of_plotted_res,dict_of_closest_res)
