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
	"DOM_FILE_4XP1",
	"AMINO_DIAGRAM",
	"DOMAIN_DIAGRAM",
	"TOPOLOGY",
	"TRAJ_20FR",
	"TRAJ_50FR",
	"CLOCK_280",
	"MOL_SVG_RMSF",
	"AMINO_RMSF_2TRAJ"
 
]

from pkg_resources import resource_filename

PDB = resource_filename(__name__,'data/4XP1.pdb')
LIG_PDB = resource_filename(__name__,'data/LIG.pdb')
MOL_SVG_INIT = resource_filename(__name__,'data/molecule_init.svg')
DOM_FILE_4XP1 = resource_filename(__name__,'data/domain_file_4XP1.txt')
AMINO_DIAGRAM = resource_filename(__name__,'data/amino_diagrams.svg')
DOMAIN_DIAGRAM = resource_filename(__name__,'data/domain_diagram.svg')
TOPOLOGY = resource_filename(__name__,'data/topology.gro')
TRAJ_20FR = resource_filename(__name__,'data/trajectory_20frames.xtc')
TRAJ_50FR = resource_filename(__name__,'data/trajectory_50frames.xtc')
CLOCK_280 = resource_filename(__name__,'data/280.svg')
MOL_SVG_RMSF = resource_filename(__name__,'data/molecule_with_rmsf.svg')
AMINO_RMSF_2TRAJ = resource_filename(__name__,'data/amino_rmsf_2traj.svg')
# This should be the last line: clean up namespace
del resource_filename
