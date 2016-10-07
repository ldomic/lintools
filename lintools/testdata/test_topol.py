from numpy.testing import TestCase, assert_equal, assert_almost_equal
import unittest
import os
import numpy as np
import MDAnalysis
from lintools.lintools import Lintools
from lintools.lintools.testdata.datafiles import *

################ AMI in P-gp with -ro 30 #################

class TestBasic1(TestCase):
    def setUp(self):
        u =MDAnalysis.Universe(AMI_GRO)
        lig_name = u.select_atoms("resname UNK")
        self.test_svg = TEST1
        self.output_name = "test1_1"
        self.lintools = Lintools(AMI_GRO,[],AMI_MOL2,lig_name,30,3.5,[1],[None],[None],0.3,"amino",self.output_name)
        self.lintools.save_files()
        self.lintools.data_input_and_res_time_analysis()
    	self.lintools.analysis_of_prot_lig_interactions()
    	self.lintools.plot_residues()
    	self.lintools.draw_figure()
    def tearDown(self):
        self.lintools.remove_files()
        del self.lintools

    def test_ami_amino(self):
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

