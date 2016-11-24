from argparse import ArgumentParser
import os
import sys
import MDAnalysis
from data import Data
from plots import Plots
from molecule import Molecule
from figure import Figure
from analysis.hbonds import HBonds
from analysis.residence_time import Residence_time
from analysis.rmsf import RMSF_measurements
from analysis.salt_bridges import SaltBridges
from analysis.pistacking import PiStacking
from analysis.sasa import SASA
from draw import Draw
from timeit import default_timer as timer
from ligand_description import LigDescr

class Lintools(object):
    """This class controls the behaviour of all other classes (Data,Plots,Molecule,Figure) 
     of lintools and inherits and transfers them resulting in a final SVG file that contains
     the protein-ligand interactions.

    It also controls the analysis (Residence_time and HBonds classes).

    Takes:
        * topology * - topology file
        * trajectory * - trajectory file(s)
        * mol2_file * - MOL2 file of the ligand
        * ligand * - MDAnalysis atomgroup of ligand that is going to be analysed
        * offset * - residue offset which determines by how many numbers the protein residue numbering
        should be offset (e.g. with offset = 30 the first residue will be changed from 1 to 30, 2 - 31, etc.)
        * cutoff * - cutoff distance in angstroms that defines the native contacts (default - 3.5A)   
        * start_frame * - start frame(s) for trajectory analysis (can be different for each trajectory)
        * end_frame * - end frame(s) for trajectory analysis (can be different for each trajectory)
        * skip * - number of frames to skip (can be different for each trajectory)
        * analysis_cutoff * - a fraction of time a residue has to fullfil the analysis parameters for (default - 0.3)
        * diagram_type * - string of the selected diagram type (e.g. "amino" or "clocks")
        * output_name * - name of the folder with results and the final SVG file

    """
    __version__ = "09.2016"
    def __init__(self,topology,trajectory,mol2_file,ligand,offset,cutoff,start_frame,end_frame,skip,analysis_cutoff,diagram_type,output_name):
        """Defines the input variables."""
        self.topology = os.path.abspath(topology)
        try:
            self.trajectory = []
            for traj in trajectory:
                self.trajectory.append(os.path.abspath(traj))
        except Exception:
            self.trajectory = []
        if mol2_file!=None:
            self.mol2_file = os.path.abspath(mol2_file)
        else:
            self.mol2_file = mol2_file
        self.ligand = ligand
        self.offset = offset
        self.cutoff = cutoff
        self.start = [None if start_frame==[None] else int(start_frame[i]) for i in range(len(trajectory))]
        self.end = [None if end_frame==[None] else int(end_frame[i]) for i in range(len(trajectory))]
        self.skip = [None if skip==[None] else int(skip[i]) for i in range(len(trajectory))]
        self.analysis_cutoff = analysis_cutoff
        self.diagram_type = diagram_type
        self.output_name = output_name
    def data_input_and_res_time_analysis(self):
        """
        Loads the data into Data() - renumbers the residues, imports mol2 file in rdkit. 
        If there are trajectories to analyse, the residues that will be plotted are determined 
        from Residence_time() analysis.
        """
        self.topol_data = Data()
        self.topol_data.load_data(self.topology,self.mol2_file,self.ligand,self.offset)
        if len(self.trajectory) == 0:
            self.topol_data.analyse_topology(self.topology,self.cutoff)
        else:
            self.res_time = Residence_time(self.topol_data,self.trajectory, self.start, self.end, self.skip,self.topology,self.mol2_file, self.ligand,self.offset)
            self.res_time.measure_residence_time(self.cutoff)
            self.res_time.define_residues_for_plotting_traj(self.analysis_cutoff)
            self.topol_data.find_the_closest_atoms(self.topology)
    def analysis_of_prot_lig_interactions(self,hydr_bonds):
        """
        The classes and function that deal with protein-ligand interaction analysis.
        """
        if hydr_bonds!=True:
            self.hbonds = HBonds(self.topol_data,self.trajectory,self.start,self.end,self.skip,self.analysis_cutoff,distance=3)
        else:
            self.hbonds=None
        self.pistacking = PiStacking(self.topol_data,self.trajectory,self.start,self.end,self.skip, self.analysis_cutoff)
        self.sasa = SASA(self.topol_data,self.trajectory)
    
        if self.trajectory!=[]:
            self.rmsf = RMSF_measurements(self.topol_data,self.topology,self.trajectory,self.ligand)
        else:
            self.rmsf=None

        self.lig_descr = LigDescr(self.topol_data,self.rmsf,self.sasa)

        self.salt_bridges = SaltBridges(self.topol_data,self.trajectory,self.lig_descr,self.start,self.end,self.skip,self.analysis_cutoff)
    def plot_residues(self):
        """
        Calls Plot() that plots the residues with the required diagram_type.
        """
        self.plots = Plots(self.topol_data,self.diagram_type)
    def draw_figure(self):
        """
        Draws molecule through Molecule() and then puts the final figure together with
        Figure().
        """
        self.molecule = Molecule(self.topol_data)
        
        self.draw = Draw(self.topol_data,self.molecule,self.hbonds,self.pistacking,self.salt_bridges,self.lig_descr)


        self.figure = Figure(self.molecule,self.topol_data,self.draw)
        self.figure.add_bigger_box()
        self.figure.manage_the_plots()
        self.figure.draw_white_circles()
        self.figure.put_everything_together()
        self.figure.write_final_draw_file(self.output_name)
    def save_files(self):
        """Saves all output from LINTools run in a single directory named after the output name."""
        os.system("mkdir "+self.output_name)
        self.workdir = os.getcwd()
        os.chdir(self.workdir+"/"+self.output_name)
    def remove_files(self):
        """Removes intermediate files."""
        file_list = ["molecule.svg","LIG.pdb"]
        for residue in self.topol_data.dict_of_plotted_res.keys():
            file_list.append(residue[1]+residue[2]+".svg")
        for f in file_list:
            if os.path.isfile(f)==True:
                os.remove(f)


if __name__ == '__main__': 
    #################################################################################################################

    parser = ArgumentParser(description='Analysis and visualisation tool for protein ligand interactions. Requires rdkit, shapely, MDAnalysis.')
    parser.add_argument('-t', '--topology', dest = 'topology', default=None, help='Topology file')
    parser.add_argument('-x', '--trajectory', dest = "trajectory", nargs="*", default=[], help='Trajectory file(s)')
    parser.add_argument('-o', '--outname', dest = "output_name", help='Name for output folder and file')
    parser.add_argument('-mol2', dest='mol2', default=None, help="MOL2 file of the ligand.")
    parser.add_argument('-c', '--cutoff', dest = "cutoff", default = 3.5, help='Cutoff distance in angstroms.')
    parser.add_argument('-ro', '--residueoffset', dest = "offset", default = 0, help='Number by which to change protein residue numbers')
    parser.add_argument('-ac', '--analysis_cutoff', dest = "analysis_cutoff", default=0.3, help='Analysis cutoff - a feature has to appear for at least a fraction of the simulation to be plotted.')
    parser.add_argument('-b', dest = "start", nargs="*", default=[None], help='Start frame number(s)')
    parser.add_argument('-e', dest = "end", nargs="*", default=[None], help='End frame number(s)')
    parser.add_argument('-skip', dest = "skip", nargs="*", default=[None], help='Skip frames')
    parser.add_argument('--no_hbonds', dest='hydr_bonds', action="store_true", help="The hydrogen bonds will not be detected.")



    args = parser.parse_args()

    ####################################################################################################################

    assert len(args.topology) >0, "No topology file provided for analysis."
    assert len(args.output_name)>0,"No output name provided."

    def find_ligand_name():
        """Users select a ligand to analyse from a numbered list."""
        gro = MDAnalysis.Universe(args.topology)
        list_of_non_ligands=["SOL","NA","CL","HOH","ARG","LYS","HIS","ASP","GLU","SER","THR", "ASN","GLN","PHE","TYR","TRP","CYS","GLY","PRO","ALA","VAL","ILE","LEU","MET"]
        potential_ligands={}
        i=0
        for residue in gro.residues:
            if residue.atoms.resnames[0] not in list_of_non_ligands:
                if residue.atoms.altLocs[0]==str("") or  residue.atoms.altLocs[0]==None:
                    potential_ligands[i]=residue.atoms
                else:
                    #Deal with ligands that have alternative locations
                    altloc = str(residue.atoms.altLocs[1])
                    resid = residue.atoms.resids[0]
                    new_residue = residue.select_atoms("resid "+str(resid)+" and altloc "+str(altloc))
                    potential_ligands[i] = new_residue
                i+=1

        print "# Nr  # Name   # Resnumber  # Chain ID"
        for lig in potential_ligands:           
            print lig, potential_ligands[lig].resnames[0], potential_ligands[lig].resids[0], potential_ligands[lig].segids[0]

        while True:
            raw = raw_input( "Choose a ligand to analyse:")
            try:
                if int(raw) in [x[0] for x in enumerate(potential_ligands.keys())] :
                    break
                else:
                    print "Error. No such group "+str(raw)
            except ValueError:
                print "Error. No such group "+str(raw)
                pass
        ligand_name=potential_ligands[int(raw)]
        return ligand_name

    def find_diagram_type():
        """User selects diagram type for the residue plots."""
        available_diagrams={1:"amino", 2:"domains",3:"clock"}
        for diagram in available_diagrams:
            print diagram, " : ", available_diagrams[diagram]
        while True:
            raw_d = raw_input( "Choose diagram type:")
            try:
                if int(raw_d)-1 in [x[0] for x in enumerate(available_diagrams.keys())] :
                    break
                else:
                    print "Error. No such group "+str(raw_d)
            except ValueError:
                print "Error. No such group "+str(raw_d)
                pass
        diagram_type=available_diagrams[int(raw_d)]
        return diagram_type
    ligand_name = find_ligand_name()
    diagram_type = find_diagram_type()
    start = timer()
    lintools = Lintools(args.topology,args.trajectory,args.mol2,ligand_name,args.offset,args.cutoff,args.start,args.end,args.skip,args.analysis_cutoff,diagram_type,args.output_name)
    lintools.save_files()
    lintools.data_input_and_res_time_analysis()
    lintools.analysis_of_prot_lig_interactions(args.hydr_bonds)
    lintools.plot_residues()
    lintools.draw_figure()
    lintools.remove_files()
    end = timer()
    print "Total time for analysis: "+str(end-start)

