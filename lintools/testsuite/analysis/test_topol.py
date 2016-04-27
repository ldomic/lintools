from numpy.testing import TestCase, assert_equal, assert_almost_equal
import unittest
import os
from lintools.lintools.testsuite.datafiles import *
import numpy as np
import MDAnalysis
from lintools.lintools.lintools import Lintools

############    Test 1 - 4xp1 with LDP amino ###############

class TestBasic1(TestCase):
    def setUp(self):
        u =MDAnalysis.Universe(FILE1)
        lig_name = u.select_atoms("resname LDP")
        self.output_name ="test1"
        self.test_svg = TEST1_SVG
        self.lintools = Lintools(FILE1,None,lig_name,0,3.5,30,"amino",None,False,False,False,"test1")
        self.lintools.get_info_about_input_and_analyse()
        self.lintools.plot_residues()
        self.lintools.draw_molecule_and_figure(tests=True)
        self.lintools.write_config_file()
    def tearDown(self):
        self.lintools.remove_files()
        del self.lintools
        file_list = ["test1.svg","test1_config.txt"]
        for f in file_list:
            if os.path.isfile(f)==True:
                os.remove(f)
    def test_4xp1_amino(self):
        # Is the molecule mol2 file produced?
        assert_equal(os.path.isfile("LIG_test.mol2"),True)
        #Is the final svg file produced?
        assert_equal(os.path.isfile(self.output_name+".svg"),True)
        if os.path.isfile(self.output_name+".svg") == True:
            with open(self.test_svg,"r") as test:
                testlines = test.readlines()
            with open(self.output_name+".svg","r") as output:
                lines = output.readlines()
                i=0
                for line in lines:
                    assert_equal(testlines[i],lines[i])
                    i+=1



class TestBasic2(TestCase):
    def setUp(self):
        u =MDAnalysis.Universe(FILE1)
        lig_name = u.select_atoms("resname LDP")
        self.output_name ="test2"
        self.test_svg = TEST2_SVG
        self.lintools = Lintools(FILE1,None,lig_name,0,3.5,30,"domains",DOM_FILE_4XP1,False,False,False,"test2")
        self.lintools.get_info_about_input_and_analyse()
        self.lintools.plot_residues()
        self.lintools.draw_molecule_and_figure(tests=True)
        self.lintools.write_config_file()
    def tearDown(self):
        self.lintools.remove_files()
        del self.lintools
        file_list = ["test2.svg","test2_config.txt"]
        for f in file_list:
            if os.path.isfile(f)==True:
                os.remove(f)
    def test_4xp1_domains(self):
        # Is the molecule mol2 file produced?
        assert_equal(os.path.isfile("LIG_test.mol2"),True)
        #Is the final svg file produced?
        assert_equal(os.path.isfile(self.output_name+".svg"),True)
        if os.path.isfile(self.output_name+".svg") == True:
            with open(self.test_svg,"r") as test:
                testlines = test.readlines()
            with open(self.output_name+".svg","r") as output:
                lines = output.readlines()
                i=0
                for line in lines:
                    assert_equal(testlines[i],lines[i])
                    i+=1

class TestBasic3(TestCase):
    def setUp(self):
        u =MDAnalysis.Universe(FILE2)
        lig_name = u.select_atoms("resname UNK")
        self.output_name ="test3"
        self.test_svg = TEST3_SVG
        self.lintools = Lintools(FILE2,[TRAJ_20_FR,TRAJ_50_FR],lig_name,30,3.5,50,"clock",None,False,False,False,"test3")
        self.lintools.get_info_about_input_and_analyse()
        self.lintools.plot_residues()
        self.lintools.draw_molecule_and_figure(tests=True)
        self.lintools.write_config_file()
    def tearDown(self):
        self.lintools.remove_files()
        del self.lintools
        file_list = ["test32.svg","test3_config.txt"]
        for f in file_list:
            if os.path.isfile(f)==True:
                os.remove(f)
    def test_two_trajectories_clock(self):
        # Is the molecule mol2 file produced?
        assert_equal(os.path.isfile("LIG_test.mol2"),True)
        #Is the final svg file produced?
        assert_equal(os.path.isfile(self.output_name+".svg"),True)
        if os.path.isfile(self.output_name+".svg") == True:
            with open(self.test_svg,"r") as test:
                testlines = test.readlines()
            with open(self.output_name+".svg","r") as output:
                lines = output.readlines()
                i=0
                for line in lines:
                    assert_equal(testlines[i],lines[i])
                    i+=1

class TestBasic4(TestCase):
    def setUp(self):
        u =MDAnalysis.Universe(FILE2)
        lig_name = u.select_atoms("resname UNK")
        self.output_name ="test4"
        self.test_svg = TEST4_SVG
        self.lintools = Lintools(FILE2,[TRAJ_20_FR,TRAJ_50_FR],lig_name,30,3.5,50,"clock",None,False,True,False,"test4")
        self.lintools.get_info_about_input_and_analyse()
        self.lintools.plot_residues()
        self.lintools.draw_molecule_and_figure(tests=True)
        self.lintools.write_config_file()
    def tearDown(self):
        self.lintools.remove_files()
        del self.lintools
        file_list = ["test34.svg","test4_config.txt"]
        for f in file_list:
            if os.path.isfile(f)==True:
                os.remove(f)
    def test_two_trajectories_clock(self):
        # Is the molecule mol2 file produced?
        assert_equal(os.path.isfile("LIG_test.mol2"),True)
        #Is the final svg file produced?
        assert_equal(os.path.isfile(self.output_name+".svg"),True)
        if os.path.isfile(self.output_name+".svg") == True:
            with open(self.test_svg,"r") as test:
                testlines = test.readlines()
            with open(self.output_name+".svg","r") as output:
                lines = output.readlines()
                i=0
                for line in lines:
                    if i == 1166: # This line contains an url which changes with every itiaration and cannot be tested
                        continue
                    else:
                        assert_equal(testlines[i],lines[i])
                    i+=1