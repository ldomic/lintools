"""
Location of data files for LINtools unit testing

Data are stored in the "data/" subdirectory
Use as :
from lintools.testsuite.datafiles import *

"""

__all__ = [
	"PDB",
	"LIG_PDB",
	"MOL2"
        "MOL_SVG_INIT"
 
]

from pkg_resources import resource_filename

PDB = resource_filename(__name__,'data/4XP1.pdb')
MOL2 = resource_filename(__name__,'data/LDP.mol2')
LIG_PDB = resource_filename(__name__,'data/LIG.pdb')
MOL_SVG_INIT = resource_filename(__name__,'data/molecule_init.svg')
# This should be the last line: clean up namespace
del resource_filename
