from argparse import ArgumentParser
import os
import shutil
import sys
import yaml
import MDAnalysis
from data import Data
from plots import Plots
from molecule import Molecule
from draw import Draw
from figure import Figure
from analysis.hbonds import HBonds
from analysis.residence_time import Residence_time
from analysis.rmsf import RMSF_measurements
from analysis.salt_bridges import SaltBridges
from analysis.pistacking import PiStacking
from analysis.sasa import SASA
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
        * mol_file * - MOL file of the ligand
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
    def __init__(self,topology,trajectory,mol_file,ligand,offset,cutoff,start_frame,end_frame,skip,analysis_cutoff,diagram_type,output_name,cfg):
        """Defines the input variables."""
        self.topology = os.path.abspath(topology)
        try:
            self.trajectory = []
            for traj in trajectory:
                self.trajectory.append(os.path.abspath(traj))
        except Exception:
            self.trajectory = []
        if mol_file!=None:
            self.mol_file = os.path.abspath(mol_file)
        else:
            self.mol_file = mol_file
        self.ligand = ligand
        self.offset = offset
        self.cutoff = cutoff
        if cfg==False:
            self.start = [None if start_frame==[None] else int(start_frame[i]) for i in range(len(trajectory))]
            self.end = [None if end_frame==[None] else int(end_frame[i]) for i in range(len(trajectory))]
            self.skip = [None if skip==[None] else int(skip[i]) for i in range(len(trajectory))]
        else:
            self.start = start_frame
            self.end = end_frame
            self.skip = skip
        self.analysis_cutoff = analysis_cutoff
        self.diagram_type = diagram_type
        self.output_name = output_name
    def data_input_and_res_time_analysis(self):
        """
        Loads the data into Data() - renumbers the residues, imports mol file in rdkit.
        If there are trajectories to analyse, the residues that will be plotted are determined
        from Residence_time() analysis.
        """
        self.topol_data = Data()
        self.topol_data.load_data(self.topology,self.mol_file,self.ligand,self.offset)
        if len(self.trajectory) == 0:
            self.topol_data.analyse_topology(self.topology,self.cutoff)
        else:
            self.res_time = Residence_time(self.topol_data,self.trajectory, self.start, self.end, self.skip,self.topology, self.ligand,self.offset)
            self.res_time.measure_residence_time(self.cutoff)
            self.res_time.define_residues_for_plotting_traj(self.analysis_cutoff)
            self.topol_data.find_the_closest_atoms(self.topology)
    def analysis_of_prot_lig_interactions(self):
        """
        The classes and function that deal with protein-ligand interaction analysis.
        """
        self.hbonds = HBonds(self.topol_data,self.trajectory,self.start,self.end,self.skip,self.analysis_cutoff,distance=3)
        self.pistacking = PiStacking(self.topol_data,self.trajectory,self.start,self.end,self.skip, self.analysis_cutoff)
        self.sasa = SASA(self.topol_data,self.trajectory)
        self.lig_descr = LigDescr(self.topol_data)
        if self.trajectory!=[]:
            self.rmsf = RMSF_measurements(self.topol_data,self.topology,self.trajectory,self.ligand,self.start,self.end,self.skip)
        self.salt_bridges = SaltBridges(self.topol_data,self.trajectory,self.lig_descr,self.start,self.end,self.skip,self.analysis_cutoff)
    def plot_residues(self):
        """
        Calls Plot() that plots the residues with the required diagram_type.
        """
        self.plots = Plots(self.topol_data,self.diagram_type)
    def draw_figure(self,data_for_color=None, data_for_size=None, data_for_clouds=None, rot_bonds=None, color_for_clouds="Blues", color_type_color="viridis"):
        """
        Draws molecule through Molecule() and then puts the final figure together with
        Figure().
        """
        self.molecule = Molecule(self.topol_data)

        self.draw = Draw(self.topol_data,self.molecule,self.hbonds,self.pistacking,self.salt_bridges,self.lig_descr)
        self.draw.draw_molecule(data_for_color, data_for_size, data_for_clouds, rot_bonds, color_for_clouds, color_type_color)

        self.figure = Figure(self.molecule,self.topol_data,self.draw)
        self.figure.add_bigger_box()
        self.figure.manage_the_plots()
        self.figure.draw_white_circles()
        self.figure.put_everything_together()
        self.figure.write_final_draw_file(self.output_name)

    def save_files(self):
        """Saves all output from LINTools run in a single directory named after the output name."""
        while True:
            try:
                os.mkdir(self.output_name)
            except Exception as e:
                self.output_name = raw_input("This directory already exists - please enter a new name:")
            else:
                break
        self.workdir = os.getcwd()
        os.chdir(self.workdir+"/"+self.output_name)
    def write_config_file(self, cfg):
        if cfg!=None:
            #copy the config file to results directory
            shutil.copy("../"+cfg, "lintools.config")
        else:
            #If there was no config file, write one
            cfg_dir = {'input':{
                            'topology':self.topology,
                            'trajectory':self.trajectory,
                            'mol file':self.mol_file,
                            'ligand':self.ligand,
                            'traj start':self.start,
                            'traj end':self.end,
                            'traj skip': self.skip,
                            'offset': self.offset,
                            'distance cutoff': self.cutoff,
                            'analysis cutoff': self.analysis_cutoff,
                            'diagram type': self.diagram_type,
                            'output name': self.output_name},
                    'representation':{
                            'data to show in color':None,
                            'data to show as size':None,
                            'data to show as cloud':None,
                            'rotatable bonds':None,
                            'cloud color scheme':'Blues',
                            'atom color scheme':'viridis',
                            'clock color scheme':'summer'}
                            }

            with open("lintools.config","wb") as ymlfile:
                yaml.dump(cfg_dir,ymlfile,default_flow_style=False)
    def remove_files(self):
        """Removes intermediate files."""
        file_list = ["molecule.svg","lig.pdb","HIS.pdb","PHE.pdb","TRP.pdb","TYR.pdb","lig.mol"]
        for residue in self.topol_data.dict_of_plotted_res.keys():
            file_list.append(residue[1]+residue[2]+".svg")
        for f in file_list:
            if os.path.isfile(f)==True:
                os.remove(f)


if __name__ == '__main__':
    #################################################################################################################

    parser = ArgumentParser(description='Analysis and visualisation tool for protein ligand interactions. Requires rdkit, shapely, MDAnalysis.')
    parser.add_argument('-cfg', '--config', dest = 'config', default=None, help='Configuration file')
    parser.add_argument('-t', '--topology', dest = 'topology', default=None, help='Topology file')
    parser.add_argument('-x', '--trajectory', dest = "trajectory", nargs="*", default=[], help='Trajectory file(s)')
    parser.add_argument('-o', '--outname', dest = "output_name", help='Name for output folder and file')
    parser.add_argument('-c', '--cutoff', dest = "cutoff", default = 3.5, help='Cutoff distance in angstroms.')
    parser.add_argument('-ac', '--analysis_cutoff', dest = "analysis_cutoff", default=0.3, help='Analysis cutoff - a feature has to appear for at least a fraction of the simulation to be plotted.')

    args = parser.parse_args()

    ####################################################################################################################



    if args.config!=None:
        #If config file exists, args.will be ingnored
        print "#####################################################################"
        print "WARNING"
        print "The arguments from command line will be ignored,"
        print "if you want to make changes, do so in the configuration file."
        print "  "
        print "  "
        print "######################################################################"
        with open(args.config, "r") as ymlfile:
            cfg = yaml.load(ymlfile)
        ## Check config file input - mainly topology and output file, also handling bad input

        lintools = Lintools(cfg['input']['topology'],cfg['input']['trajectory'],cfg['input']['mol file'],cfg['input']['ligand'],cfg['input']['offset'],float(cfg['input']['distance cutoff']),cfg['input']['traj start'],cfg['input']['traj end'],cfg['input']['traj skip'],cfg['input']['analysis cutoff'],cfg['input']['diagram type'],cfg['input']['output name'],cfg=True)
        lintools.save_files()
        lintools.data_input_and_res_time_analysis()
        lintools.analysis_of_prot_lig_interactions()
        lintools.plot_residues(cfg['representation']['clock color scheme'])
        lintools.write_config_file(args.config)
        lintools.draw_figure(cfg['representation']['data to show in color'], cfg['representation']['data to show as size'], cfg['representation']['data to show as cloud'], cfg['representation']['rotatable bonds'],cfg['representation']['cloud color scheme'], cfg['representation']['atom color scheme'])
        lintools.remove_files()
    else:
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
                    try:
                        if residue.atoms.altLocs[0]==str("") or  residue.atoms.altLocs[0]==None:
                            potential_ligands[i]=residue.atoms
                        else:
                            #Deal with ligands that have alternative locations
                            altloc = str(residue.atoms.altLocs[1])
                            resid = residue.atoms.resids[0]
                            new_residue = residue.select_atoms("resid "+str(resid)+" and altloc "+str(altloc))
                            potential_ligands[i] = new_residue
                    except Exception as e:
                        potential_ligands[i]=residue.atoms
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
            return "resid "+str(ligand_name.resids[0])+" and segid "+str(ligand_name.segids[0])

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

        lintools = Lintools(args.topology,args.trajectory,None,ligand_name,0,args.cutoff,[None],[None],[None],float(args.analysis_cutoff),diagram_type,args.output_name,cfg=False)
        lintools.save_files()
        lintools.data_input_and_res_time_analysis()
        lintools.analysis_of_prot_lig_interactions()
        lintools.plot_residues()
        lintools.write_config_file(None)
        lintools.draw_figure()
        lintools.remove_files()
