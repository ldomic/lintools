from numpy.testing import TestCase, assert_equal, assert_almost_equal
import unittest
import os
from lintools.topol import Topol_Data
from lintools.testsuite.datafiles import *
import numpy as np
from lintools.molecule import Molecule

class TestCheckMolecule(TestCase):
    "Tests without running Hydrogen bond analysis as that changes the initial setup"
    def setUp(self):
        self.topology = Topol_Data(PDB,None,None,0,MOL2)
        self.u = self.topology.universe
        self.topology.ligand = self.u.select_atoms("resname LDP")
        self.topology.ligand_no_H=self.u.select_atoms("resname LDP and not name H*")
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
        init_res_plac = {'ALA117': 1479.1775886844346,
 'ASP121': -86.62591647218265,
 'ASP46': -1348.3439000483465,
 'PHE325': -343.2044430298675,
 'PHE43': -717.1381088829888,
 'SER421': -1509.0302300424964,
 'SER422': 50.009642756242464,
 'TYR124': -811.9622777596875,
 'VAL120': -348.318054709263}
        self.molecule.load_molecule_in_rdkit_smiles()
        self.molecule.convex_hull()
        assert_equal(init_res_plac,self.molecule.b_for_all)
    def test_multiple_hulls(self):
        self.molecule.load_molecule_in_rdkit_smiles()
        self.molecule.convex_hull()
        self.molecule.make_multiple_hulls()
        nearest_point_coords = {'ALA117': (496.3834752784929, 392.4554597280188),
 'ASP121': (659.5310248283491, 96.0861185646005),
 'ASP46': (-408.9438685564142, 66.88155992276614),
 'PHE325': (728.5595499673886, 349.2156810524641),
 'PHE43': (97.33123250090364, 363.80648002850836),
 'SER421': (715.8712730749303, 161.08185674420045),
 'SER422': (476.6654077739911, -107.18624472494676),
 'TYR124': (30.724337778050426, 310.7951910294776),
 'VAL120': (678.2026350259333, 354.5832074723733)}
        assert_equal(nearest_point_coords,self.molecule.nearest_points_coords)
    def test_placement(self):
        "Will have to change after update."
        self.molecule.load_molecule_in_rdkit_smiles()
        self.molecule.convex_hull()
        self.molecule.make_multiple_hulls()
        self.molecule.make_new_projection_values()
        nearest_point_coords = {'ALA117': (666.4425252891505, 180.69220611610368),
 'ASP121': (563.4877978766755, -15.018524395081442),
 'ASP46': (38.10666676729852, 3.1652860579481654),
 'PHE325': (684.3682747070243, 343.3230590400053),
 'PHE43': (72.13746038161935, 348.5459952862316),
 'SER421': (729.3476730705256, 244.81103798217305),
 'SER422': (480.86830509649735, -89.15153775487737),
 'TYR124': (-3.3678074246929555, 290.30903949338517),
 'VAL120': (598.7941318831142, 371.1989308284794)}

        assert_equal(self.molecule.nearest_points_coords, nearest_point_coords)
        

