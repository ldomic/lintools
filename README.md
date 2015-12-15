# LINTools

LINTools (or Ligand Interaction Network Tools) is an open source program fully written in Python. It produces diagrams of interaction between protein and ligands in a molecular dynamics trajectory or a topology file. It can be used to investigate docking poses, show the residues that are spending most time in the vicinity of the ligand and possibly other things. Examples shall follow on a later date.

#Dependencies
LINTools require these packages:
* RDKit (instalation instructions with Anaconda available from RDKit GitHub account: https://github.com/rdkit/conda-rdkit)
* Shapely (available on GitHub https://github.com/Toblerity/Shapely)
* MDAnalysis (available on GitHub https://github.com/MDAnalysis/mdanalysis)

If this is a problem for your computer's architecture, a Dockerfile has also been provided.

You will need at least a topology file and a mol2 file of your ligand.

Usage:
For a topology file (no trajectory data):
```
python lintools.py -t my_top_file.pdb -o my_output -m ligand.mol2
(Optional: --cutoff [a number] --residueoffset [a number] --diagram_type "amino" )
```

For trajectory data:
```
python lintools.py -t my_top_file.pdb -x my_traj.xtc -o my_output -m ligand.mol2
(Optional: same as above + --analysis "occurance"or  "rmsf" and for occurance analysis it is possible to choose up to three trajectories
which are displayed as clock diagrams. --diagram_type "clock" must be specified)
```
