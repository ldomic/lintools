from numpy.testing import TestCase, assert_equal, assert_almost_equal
import unittest
import os
from lintools.lintools.topol import Topol_Data
from lintools.lintools.plots import Plots
from lintools.lintools.molecule import Molecule
from lintools.lintools.analysis.hbonds import HBonds
from lintools.lintools.figure import Figure
from lintools.lintools.testsuite.datafiles import *
import numpy as np

class TestCheckFigure(TestCase):
    def setUp(self):
        self.topology = Topol_Data(PDB,None,None,0)
        self.u = self.topology.universe
        self.ligand = self.u.select_atoms("resname LDP")
        self.ligand.resname = "LIG"
        self.ligand.resnames = "LIG"
        self.topology.define_ligand(self.ligand)
        self.topology.find_res_to_plot(3.5)
        self.hbonds = HBonds(self.topology,PDB,None,self.ligand,0,30)
        self.topology.get_closest_ligand_atoms(self.hbonds)
        self.molecule = Molecule(self.topology)
        self.plots = Plots(self.topology)
    def tearDown(self):
    	del self.topology
    	del self.hbonds
    	del self.molecule
    	del self.plots
    	del self.figure
    	if os.path.isfile("amino_diagrams.svg")==True:
            os.remove("amino_diagrams.svg")
    def test_plot_amino_diagrams(self):
    	self.plots.define_amino_acids()
    	self.plots.plot_amino_diagramms()
    	self.figure = Figure(self.molecule,"amino",self.topology,self.hbonds,self.plots,tests=True)
        self.figure.draw_hbonds_in_graph()
        self.figure.draw_white_circles_at_atoms()
        self.figure.put_everything_together()
        self.figure.write_final_draw_file("amino_diagrams")
        with open(AMINO_DIAGRAM,"r") as f:
            lines = f.readlines()
            self.out_test_svg = " ".join(map(str,lines[2:-1]))
            f.close()
        with open("amino_diagrams.svg","r") as f:
            lines = f.readlines()
            self.out_svg_to_test = " ".join(map(str,lines[2:-1]))
            f.close()
        assert_equal(self.out_test_svg,self.out_svg_to_test)
    def test_plot_domains_diagrams(self):
        self.plots.define_domains(DOM_FILE_4XP1, 3.5)
        self.plots.plot_domains_diagramms()
        self.figure = Figure(self.molecule,"domain",self.topology,self.hbonds,self.plots,tests=True)
        self.figure.draw_hbonds_in_graph()
        self.figure.draw_white_circles_at_atoms()
        self.figure.put_everything_together()
        self.figure.write_final_draw_file("domain_diagrams")
        with open(DOMAIN_DIAGRAM,"r") as f:
            lines = f.readlines()
            self.out_test_svg = " ".join(map(str,lines[2:-1]))
            f.close()
        with open("domain_diagrams.svg","r") as f:
            lines = f.readlines()
            self.out_svg_to_test = " ".join(map(str,lines[2:-1]))
            f.close()
        assert_equal(self.out_test_svg,self.out_svg_to_test)

