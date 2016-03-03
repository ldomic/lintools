if __name__ == '__main__': 
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
	from topol import Config

	#################################################################################################################

	parser = ArgumentParser(description='Analysis and visualisation tool for protein ligand interactions. Requires rdkit, shapely, MDAnalysis modules.')
	parser.add_argument('-t', '--topology', dest = 'grofile', default=None, help='Input File name of topology file. Accepts gro, pdb files')
	parser.add_argument('-x', '--trajectory', dest = "xtcfile", nargs="*", default=None, help='Input File name of trajectory file(s). Accepts up to 3 xtc files (Optional. Default: None)')
	parser.add_argument('-o', '--outname', dest = "output_name", help='Name of the output files.')
	parser.add_argument('-a', '--analysis', dest = "analysis_type", default = None, help='Select type of analysis for plotting. Available types - RMSF, occurrence. (Optional, default is None.)')
	parser.add_argument('-c', '--cutoff', dest = "cutoff", default = 3.5, help='Input cutoff distance from the ligand that is taken into account in angstroms (Example: 3.5).')
	parser.add_argument('-ro', '--residueoffset', dest = "offset", default = 0, help='Input the number of offset residues for the protein. (Optional, default is 0)')
	parser.add_argument('-ac', '--analysis_cutoff', dest = "analysis_cutoff", default=30, help='Analysis cutoff - a feature has to appear for at least a third of the simulation to be counted. Default: 30')
	parser.add_argument('-df', '--domain_file', dest = "domain_file", default=None, help='Input file for domains of your protein. To see the required format, check README or our GitHub page')
	parser.add_argument('-conf', '--config_file', dest = "config_file", default=None, help="Input the name of the config file.")
	parser.add_argument('--no_hbonds', dest='hydr_bonds', action="store_true", help="The hydrogen bonds will not be detected.")
	parser.add_argument('--debug', dest='debug', action="store_true", help="Functions for debugging.")


	args = parser.parse_args()

	###################################################################################################################
	
	### All assertions

	assert len(args.output_name)>0, "Provide a name for output file"
	####################################################################################################################
	###  Define input specifications:
	if args.config_file!=None:
		config_read=Config()
		config_read.read_config_file(args.config_file)
		# Topology
		if args.grofile!=None:
			topology = os.path.abspath(args.grofile)
		else:
			topology=config_read.topology
		# Trajectory
		if args.xtcfile!=None:
			i=0
			trajectory=[]
			for i in range(len(args.xtcfile)):
				trajectory.append(os.path.abspath(args.xtcfile[i]))
		else:
			if config_read.trajectory!=None:
				trajectory=config_read.trajectory
			else:
				trajectory=None
		# Offset
		if args.offset==0:
			offset = config_read.res_offset
		else:
			offset=args.offset
		# Domain file
		if args.domain_file!=None:
			domain_file=os.path.abspath(args.domain_file)
		else:
			domain_file=config_read.domain_file
		# Analysis type
		if args.analysis_type!=None:
			analysis_type= args.analysis_type
		else:
			analysis_type=config_read.analysis_type
		# Cutoff
		if args.cutoff!=3.5:
			cutoff=args.cutoff
		else:
			cutoff=config_read.cutoff
		#Analysis cutoff
		if args.analysis_cutoff!=30:
			analysis_cutoff=args.analysis_cutoff
		else:
			analysis_cutoff=config_read.analysis_cutoff
	else:
		assert len(args.grofile)>0, "Provide an input topology file"
		topology=os.path.abspath(args.grofile)
		if args.xtcfile!=None:
			i=0
			trajectory=[]
			for i in range(len(args.xtcfile)):
				trajectory.append(os.path.abspath(args.xtcfile[i]))
		else:
			trajectory=None
		offset=args.offset
		if args.domain_file!=None:
			domain_file=os.path.abspath(args.domain_file)
		else:
			domain_file=None
		#diagram_type=args.diagram_type
		analysis_type=args.analysis_type
		cutoff=args.cutoff
		analysis_cutoff=args.analysis_cutoff
 	

 	#######################################################################################################################

	gro = MDAnalysis.Universe(topology)
	list_of_non_ligands=["SOL","NA","CL","HOH","ARG","LYS","HIS","ASP","GLU","SER","THR", "ASN","GLN","PHE","TYR","TRP","CYS","GLY","PRO","ALA","VAL","ILE","LEU","MET"]
	if args.config_file!=None:
		potential_ligands={0:"Ligand molecule from config file"}
		i=1
	else:
		potential_ligands={}
		i=0
	for residue in gro.residues:
	    if residue.resnames[0] not in list_of_non_ligands:
	        potential_ligands[i]=residue
	        i+=1

	print "# Nr  # Name   # Resnumber  # Chain ID"
	if args.config_file!=None:
		for lig in potential_ligands:
			if lig==0:
				print lig, potential_ligands[lig]
			else:
				print lig, potential_ligands[lig].resnames[0], potential_ligands[lig].resids[0], potential_ligands[lig].segids[0]
	else:
		for lig in potential_ligands:			
			print lig, potential_ligands[lig].resnames[0], potential_ligands[lig].resids[0], potential_ligands[lig].segids[0]
	ligand_name=potential_ligands[int(raw_input( "Choose a ligand to analyse:"))]
	print ligand_name, ligand_name.resids[0], ligand_name.segids[0]
	if args.config_file!=None:
		if ligand_name==potential_ligands[0]:
			for ligand in potential_ligands.keys()[1:]:
				if potential_ligands[ligand].resnames[0]==config_read.ligand_name[0] and potential_ligands[ligand].resids[0]==int(config_read.ligand_name[1]) and potential_ligands[ligand].segids[0]==config_read.ligand_name[2]:
					ligand_name=potential_ligands[ligand]


	if args.config_file!=None:
		if args.domain_file!=None:
			if args.analysis_type=="occurrence" or trajectory!=None:
				available_diagrams={0:"From config file",1:"amino", 2:"domains",3:"clock"}
			if trajectory==None:
				available_diagrams={0:"From config file", 1:"amino", 2:"domains"}
		else:
			if args.analysis_type=="occurrence" or trajectory!=None:
				available_diagrams={0:"From config file",1:"amino", 2:"clock"}
			if trajectory==None:
				available_diagrams={0:"From config file",1:"amino", 2:"clock"}
		for diagram in available_diagrams:
			print diagram, " : ", available_diagrams[diagram]
		diagram_type=available_diagrams[int(raw_input("Choose diagram type:"))]
		if diagram_type=="From config file":
			diagram_type=config_read.diagram_type
	else:
		if args.domain_file!=None:
			if args.analysis_type=="occurrence" or trajectory!=None:
				available_diagrams={1:"amino", 2:"domains",3:"clock"}
			if trajectory==None:
				available_diagrams={1:"amino", 2:"domains"}
		else:
			if args.analysis_type=="occurrence" or trajectory!=None:
			#if args.analysis_type=="occurrence":
				available_diagrams={1:"amino", 2:"clock"}
			if trajectory==None:
				available_diagrams={1:"amino"}
		for diagram in available_diagrams:
			print diagram, " : ", available_diagrams[diagram]
		diagram_type=available_diagrams[int(raw_input("Choose diagram type:"))]


	#############################################################################################################

	if analysis_type=="occurrence":
		md_sim = Topol_Data(topology, None, ligand_name, offset)
		md_sim.define_ligand(ligand_name)
		md_sim.make_mol2_file()
		occurrence = Occurrence_analysis(topology, trajectory, ligand_name, cutoff, offset, md_sim)
		occurrence.get_closest_residues(analysis_cutoff)
		if args.hydr_bonds!=True:
			hbonds = HBonds(md_sim,topology, trajectory, ligand_name, offset,analysis_cutoff)
	else: 
	#if analysis type is anything different only one traj at time is going to be analysed
		assert trajectory is None or len(trajectory)<=1, "Only one trajectory at the time can be analysed."
		if trajectory	is None:
			md_sim = Topol_Data(topology, trajectory, ligand_name, offset)
			md_sim.define_ligand(ligand_name)
			md_sim.make_mol2_file()
			md_sim.find_res_to_plot(cutoff)
			if args.hydr_bonds!=True:
				hbonds = HBonds(md_sim,topology, trajectory, ligand_name, offset,analysis_cutoff)
		else:
			md_sim = Topol_Data(topology, None, ligand_name, offset)
			md_sim.define_ligand(ligand_name)
			md_sim.make_mol2_file()
			occurrence = Occurrence_analysis(topology, trajectory, ligand_name, cutoff, offset, md_sim)
			occurrence.get_closest_residues(analysis_cutoff)
			if args.hydr_bonds!=True:
				hbonds = HBonds(md_sim, topology, trajectory, ligand_name, offset,analysis_cutoff)

	if args.hydr_bonds!=True:
		md_sim.get_closest_ligand_atoms(hbonds)
	else:
		md_sim.get_closest_ligand_atoms()


	plots = Plots(md_sim)
	if diagram_type=="amino":
		plots.define_amino_acids()
		plots.plot_amino_diagramms()
	if diagram_type=="domains":
		assert len(domain_file)>0, "Provide a file defining domains"
		plots.define_domains(domain_file, offset)
		plots.plot_domains_diagramms()
	if diagram_type=="clock":
		plots.plot_clock_diagramms()


	molecule = Molecule(md_sim)
	if analysis_type=="RMSF" or analysis_type=="rmsf":
		rmsf = RMSF_measurements(md_sim,topology, trajectory, ligand_name, offset, args.output_name)
		molecule = Molecule(md_sim, rmsf)



	if analysis_type=="RMSF" or analysis_type=="rmsf":
		if args.hydr_bonds!=True:
			figure=Figure(molecule, diagram_type,md_sim,hbonds,plots,rmsf)
		else:
			figure=Figure(molecule, diagram_type,md_sim, plot_object=plots,rmsf_object=rmsf)
	else:
		if args.hydr_bonds!=True:
			figure=Figure(molecule, diagram_type,md_sim,hbonds,plots)
		else:
			figure=Figure(molecule, diagram_type,md_sim,plots)
	if args.hydr_bonds!=True:
		figure.draw_hbonds_in_graph()
	figure.draw_white_circles_at_atoms()
	if args.debug==True:
		figure.draw_lines_in_graph() #a function for debugging purposes
	figure.put_everything_together()
	figure.write_final_draw_file(args.output_name)

	file_list=[]
	file_list=["molecule.svg"]
	for residue in md_sim.dict_of_plotted_res.keys():
		file_list.append(str(residue[3:])+".svg")

	for f in file_list:
		os.remove(f)

	config_write = Config(md_sim)
	config_write.write_config_file(args.output_name, topology, trajectory, offset, diagram_type, cutoff, analysis_type, domain_file,analysis_cutoff)
	print "Ready!"

	print "The figure you created has been saved under filename "+args.output_name+".svg"
	print "This file format can be converted to any other image file format using "
	print "free software such as ImageMagick, Inkscape or GIMP."
	print "    "
	print "#########################################################################"
	print "       "
	print "Thank you for using Ligand Interaction Network. Check our GitHub account"
	print "https://github.com/ldomic/lintools."
	
