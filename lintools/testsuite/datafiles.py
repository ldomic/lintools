"""
Location of data files for LINtools unit testing

Data are stored in the "data/" subdirectory
Use as :
from lintools.testsuite.datafiles import *

"""

__all__ = [
	"PDB",
	"LIG_PDB",
        "MOL_SVG_INIT",
	"DOM_FILE_4XP1"
 
]

from pkg_resources import resource_filename

PDB = resource_filename(__name__,'data/4XP1.pdb')
LIG_PDB = resource_filename(__name__,'data/LIG.pdb')
MOL_SVG_INIT = resource_filename(__name__,'data/molecule_init.svg')
DOM_FILE_4XP1 = resource_filename(__name__,'data/domain_file_4XP1.txt')
# This should be the last line: clean up namespace
del resource_filename
