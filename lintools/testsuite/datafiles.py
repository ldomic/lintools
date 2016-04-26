"""
Location of data files for LINtools unit testing

Data are stored in the "data/" subdirectory
Use as :
from lintools.testsuite.datafiles import *

"""

__all__ = [
	"FILE1",
	"FILE2",
	"TRAJ_20_FR",
	"TRAJ_50_FR",
	"TEST1_SVG",
	"DOM_FILE_4XP1",
	"TEST2_SVG",
	"TEST3_SVG",
	"TEST4_SVG"
]

from pkg_resources import resource_filename

FILE1 = resource_filename(__name__,'data/4XP1.pdb')
DOM_FILE_4XP1 = resource_filename(__name__,'data/domain_file_4XP1.txt')
TEST1_SVG = resource_filename(__name__,'data/amino_diagram_test1.svg')
TEST2_SVG = resource_filename(__name__,'data/domain_diagram_test2.svg')
FILE2 = resource_filename(__name__,'data/topology.gro')
TRAJ_20_FR = resource_filename(__name__,'data/trajectory_20frames.xtc')
TRAJ_50_FR = resource_filename(__name__,'data/trajectory_50frames.xtc')
TEST3_SVG = resource_filename(__name__,'data/clock_diagram_test3.svg')
TEST4_SVG = resource_filename(__name__,'data/clock_rmsf_test4.svg')
# This should be the last line: clean up namespace
del resource_filename
