__version__ = '0.1-dev04.2016'
from argparse import ArgumentParser
import os
import sys
import MDAnalysis
from topol import Topol_Data
from plots import Plots
from molecule import Molecule
from figure import Figure
from analysis.hbonds import HBonds
from analysis.rmsf import RMSF_measurements
from analysis.occurrence import Occurrence_analysis

class Lintools(object):
	def __init__(self, topology,trajectory,ligand_name,offset,cutoff,analysis_cutoff,diagram_type,domain_file,HB_flag,RMSF_flag,debug_flag,output_name):
		self.topology = topology
		self.trajectory = trajectory
		self.ligand_name  = ligand_name
		self.offset = offset
		self.cutoff = cutoff
		self.analysis_cutoff = analysis_cutoff
		self.diagram_type = diagram_type
		self.domain_file = domain_file
		self.HB_flag = HB_flag
		self.RMSF_flag = RMSF_flag
		self.debug_flag = debug_flag
		self.output_name = output_name
		self.rmsf = None
		self.hbonds = None
	def get_info_about_input_and_analyse(self):
		"""This function loads all input files and decides which residues to plot"""
		self.topol_data = Topol_Data(self.topology, self.trajectory, self.ligand_name, self.offset)
		self.topol_data.define_ligand(self.ligand_name)
		if self.trajectory==None:
			self.topol_data.find_res_to_plot(self.cutoff)
		else:
			occurrence = Occurrence_analysis(self.topology, self.trajectory, self.ligand_name, self.cutoff, self.offset, self.topol_data)
			occurrence.get_closest_residues(self.analysis_cutoff)
		if self.HB_flag!=True:
			self.hbonds = HBonds(self.topol_data, self.topology, self.trajectory, self.ligand_name, self.offset,self.analysis_cutoff)
			self.topol_data.get_closest_ligand_atoms(self.hbonds)
		else:
			self.topol_data.get_closest_ligand_atoms()
		if self.RMSF_flag==True:
			self.rmsf = RMSF_measurements(self.topol_data,self.topology, self.trajectory, self.ligand_name, self.offset, self.output_name)

	def plot_residues(self):
	    self.plots = Plots(self.topol_data)
            if self.diagram_type=="amino":
		self.plots.define_amino_acids()
		self.plots.plot_amino_diagramms()
	    if self.diagram_type=="domains":
		assert len(self.domain_file)>0, "Provide a file defining domains"
		self.plots.define_domains(self.domain_file, self.offset)
		self.plots.plot_domains_diagramms()
	    if self.diagram_type=="clock":
		self.plots.plot_clock_diagramms()

	def draw_molecule_and_figure(self,tests=False):
	    self.molecule = Molecule(self.topol_data, self.rmsf)
            self.figure=Figure(self.molecule, self.diagram_type,self.topol_data,self.hbonds,self.plots,self.rmsf, tests)
            self.figure.draw_hbonds_in_graph()
            self.figure.draw_white_circles_at_atoms()
            if self.debug_flag==True:
		self.figure.draw_lines_in_graph() #a function for debugging purposes
            self.figure.put_everything_together()
            self.figure.write_final_draw_file(self.output_name)


	def remove_files(self):
            file_list = ["molecule.svg","LIG.pdb","LIG_test.mol2","test.xtc","rmsf_colorbar.svg"]
            for residue in self.topol_data.dict_of_plotted_res.keys():
                file_list.append(str(residue[3:])+".svg")
                for f in file_list:
                    if os.path.isfile(f)==True:
                        os.remove(f)


if __name__ == '__main__': 
	#################################################################################################################

	parser = ArgumentParser(description='Analysis and visualisation tool for protein ligand interactions. Requires rdkit, shapely, MDAnalysis modules.')
	parser.add_argument('-t', '--topology', dest = 'topology', default=None, help='Input File name of topology file. Accepts gro, pdb files')
	parser.add_argument('-x', '--trajectory', dest = "trajectory", nargs="*", default=None, help='Input File name of trajectory file(s). Accepts up to 3 xtc files (Optional. Default: None)')
	parser.add_argument('-o', '--outname', dest = "output_name", help='Name of the output files.')
	parser.add_argument('-rmsf', '--rmsf', dest = "rmsf", action="store_true", help="Analysis of ligand root mean square fluctuations.")
	parser.add_argument('-c', '--cutoff', dest = "cutoff", default = 3.5, help='Input cutoff distance from the ligand that is taken into account in angstroms (Example: 3.5).')
	parser.add_argument('-ro', '--residueoffset', dest = "offset", default = 0, help='Input the number of offset residues for the protein. (Optional, default is 0)')
	parser.add_argument('-ac', '--analysis_cutoff', dest = "analysis_cutoff", default=30, help='Analysis cutoff - a feature has to appear for at least a third of the simulation to be counted. Default: 30')
	parser.add_argument('-df', '--domain_file', dest = "domain_file", default=None, help='Input file for domains of your protein. To see the required format, check README or our GitHub page')
	parser.add_argument('--no_hbonds', dest='hydr_bonds', action="store_true", help="The hydrogen bonds will not be detected.")
	parser.add_argument('--debug', dest='debug', action="store_true", help="Functions for debugging.")


	args = parser.parse_args()

    ####################################################################################################################
	def find_ligand_name():
		gro = MDAnalysis.Universe(args.topology)
		list_of_non_ligands=["SOL","NA","CL","HOH","ARG","LYS","HIS","ASP","GLU","SER","THR", "ASN","GLN","PHE","TYR","TRP","CYS","GLY","PRO","ALA","VAL","ILE","LEU","MET"]
		potential_ligands={}
		i=0
		for residue in gro.residues:
		    if residue.resnames[0] not in list_of_non_ligands:
		    	if residue.altLocs[0]==str("") or  residue.altLocs[0]==None:
		        	potential_ligands[i]=residue.atoms
		    	else:
		    		#Deal with ligands that have alternative locations
		    		altloc = str(residue.altLocs[1])
		    		resid = residue.resids[0]
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
		if args.domain_file!=None:
			if  args.trajectory!=None:
				available_diagrams={1:"amino", 2:"domains",3:"clock"}
			if args.trajectory==None:
				available_diagrams={1:"amino", 2:"domains"}
		else:
			if  args.trajectory!=None:
				available_diagrams={1:"amino", 2:"clock"}
			if args.trajectory==None:
				available_diagrams={1:"amino"}
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

			###################################################################################################################
	if type(args.output_name)!=str:
		raise IOError,"Provide a name for output file"
	ligand_name = find_ligand_name()
	diagram_type = find_diagram_type()
	lintools = Lintools(args.topology, args.trajectory, ligand_name, args.offset, args.cutoff, args.analysis_cutoff, diagram_type, args.domain_file, args.hydr_bonds, args.rmsf, args.debug, args.output_name)
	lintools.get_info_about_input_and_analyse()
	lintools.plot_residues()
	lintools.draw_molecule_and_figure()
	lintools.remove_files()
