"""
Location of data files for LINtools unit testing

Data are stored in the "data/" subdirectory
Use as :
from lintools.testsuite.datafiles import *

"""

__all__ = [
	"AMI_GRO",
	"AMI_MOL2",
	"AMI_XTC",
	"TEST1",
	"PDB_4XP1",
	"LDP_MOL2",
	"TEST2",
	"TEST3"
]

from pkg_resources import resource_filename

AMI_GRO = resource_filename(__name__,'ami/lig.gro')
AMI_MOL2 = resource_filename(__name__,'ami/lig.mol2')
AMI_XTC = resource_filename(__name__,'ami/lig.xtc')
TEST1 = resource_filename(__name__,'test1/test1.svg')
PDB_4XP1 = resource_filename(__name__,'ldp/4XP1.pdb')
LDP_MOL2 = resource_filename(__name__,'ldp/LDP.mol2')
TEST2 = resource_filename(__name__,'test2/test2.svg')
TEST3 = resource_filename(__name__,'test3/test3.svg')
# This should be the last line: clean up namespace
del resource_filename