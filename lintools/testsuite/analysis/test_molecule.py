from numpy.testing import TestCase, assert_equal, assert_almost_equal
import unittest
import os
from lintools.lintools.topol import Topol_Data
from lintools.lintools.testsuite.datafiles import *
import numpy as np
from lintools.lintools.molecule import Molecule

class TestCheckMolecule(TestCase):
    "Tests without running Hydrogen bond analysis as that changes the initial setup"
    def setUp(self):
        self.topology = Topol_Data(PDB,None,None,0)
        self.u = self.topology.universe
        self.topology.ligand = self.u.select_atoms("resname LDP")
        self.topology.ligand_no_H=self.u.select_atoms("resname LDP and not name H*")
        self.topology.define_ligand(self.topology.ligand)
        self.topology.find_res_to_plot(3.5)
        self.topology.get_closest_ligand_atoms()
        self.molecule = Molecule(self.topology,test=True)
    def tearDown(self):
        del self.topology
        del self.molecule
    def test_smiles(self):
        "SMILES"
        self.molecule.load_molecule_in_rdkit_smiles()
        assert_equal = (Chem.MolToSmiles(self.molecule.smiles),"NCCc1ccc(O)c(O)c1")
    def test_molecule_svg(self):
        "Drawing of initial molecule"
        self.molecule.load_molecule_in_rdkit_smiles()
        with open(MOL_SVG_INIT,"r") as f:
            lines = f.readlines()
            self.out_test_svg = " ".join(map(str,lines[2:-1]))
            f.close()
        with open("molecule.svg","r") as f:
            lines = f.readlines()
            self.out_svg_to_test = " ".join(map(str,lines[2:-1]))
            f.close()
        assert_equal(self.out_test_svg,self.out_svg_to_test)
    def test_atom_identities(self):
        "Test whether atom identities between SMILES and PDB formats are assigned correctly."
        self.molecule.load_molecule_in_rdkit_smiles()
        atom_identities  = {'0': 2,
 '1': 3,
 '10': 0,
 '2': 6,
 '3': 10,
 '4': 4,
 '5': 5,
 '6': 8,
 '7': 9,
 '8': 7,
 '9': 1}
        assert_equal(atom_identities,self.molecule.atom_identities)
    def test_init_residue_placement(self):
        init_res_plac = {'ALA117': 2308.1994762176714,
 'ASP121': 1831.0082742001227,
 'ASP46': 472.94554763298447,
 'PHE325': 31.0582201951117,
 'PHE43': 267.5310150287892,
 'SER421': 1466.0559780001017,
 'SER422': 1831.0082742001227,
 'TYR124': 1233.1374544521545,
 'VAL120': 2410.90671796573}
        self.molecule.load_molecule_in_rdkit_smiles()
        self.molecule.convex_hull()
        assert_equal(init_res_plac,self.molecule.b_for_all)
    def test_multiple_hulls(self):
        self.molecule.load_molecule_in_rdkit_smiles()
        self.molecule.convex_hull()
        self.molecule.make_multiple_hulls()
        nearest_point_coords = {'ALA117': (905.3165439777306, 273.55858578116874),
 'ASP121': (614.122000806857, 520.449030803404),
 'ASP46': (63.02200302968498, 96.6579216550051),
 'PHE325': (472.75358992150205, -96.59143965915837),
 'PHE43': (227.40714798240575, -36.685488837522726),
 'SER421': (409.65009360496936, 530.1469317852273),
 'SER422': (604.3330816038679, 555.0271159857422),
 'TYR124': (196.95331974896874, 426.38242249470085),
 'VAL120': (806.0821639624614, -13.885096431484556)}
        assert_equal(nearest_point_coords,self.molecule.nearest_points_coords)
    def test_placement(self):
        self.molecule.load_molecule_in_rdkit_smiles()
        self.molecule.convex_hull()
        self.molecule.make_multiple_hulls()
        self.molecule.make_new_projection_values()
        nearest_point_coords = {'ALA117': (905.3165439777306, 273.55858578116874),
 'ASP121': (574.5185749015425, 509.2374575105398),
 'ASP46': (63.02200302968498, 96.6579216550051),
 'PHE325': (472.75358992150205, -96.59143965915837),
 'PHE43': (227.40714798240575, -36.685488837522726),
 'SER421': (409.65009360496936, 530.1469317852273),
 'SER422': (643.9365075091823, 566.2386892786062),
 'TYR124': (196.95331974896874, 426.38242249470085),
 'VAL120': (806.0821639624614, -13.885096431484556)}

        assert_equal(self.molecule.nearest_points_coords, nearest_point_coords)
        

