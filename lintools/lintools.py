if __name__ == '__main__': 
	from argparse import ArgumentParser
	import os
	import MDAnalysis
	from topol import Topol_Data
	from plots import Plots
	from molecule import Molecule
	from figure import Figure
	from analysis.hbonds import HBonds
	from analysis.rmsf import RMSF_measurements
	from analysis.occurance import Occurance_analysis

	#################################################################################################################

	parser = ArgumentParser(description='Analysis and visualisation tool for protein ligand interactions. Requires rdkit, shapely, MDAnalysis modules.')
	parser.add_argument('-t', '--topology', dest = 'grofile', help='Input File name of topology file. Accepts gro, pdb files')
	parser.add_argument('-x', '--trajectory', dest = "xtcfile", nargs="*", default=None, help='Input File name of trajectory file(s). Accepts up to 3 xtc files (Optional. Default: None)')
	parser.add_argument('-o', '--outname', dest = "output_name", help='Name of the output svg file.')
	parser.add_argument('-a', '--analysis', dest = "analysis_type", default = None, help='Select type of analysis for plotting. Available types - RMSF of ligand, occurance_analysis. (Optional, default is None.)')
	parser.add_argument('-c', '--cutoff', dest = "cutoff", default = 3.5, help='Input cutoff distance from the ligand that is taken into account in angstroms (Example: 3.5).')
	parser.add_argument('-ro', '--residueoffset', dest = "offset", default = 0, help='Input the number of offset residues for the protein. (Optional, default is 0)')
	parser.add_argument('-d', '--diagram_type', dest = "diagram_type", default="amino", help='Input type of diagramm required. Options: "clock" for clock diagrams (Only available with trajectory present), "helices" for diagrams representing residue membership to certain helices (requires user input to determine helices), "amino" showing the amino acid type')
	parser.add_argument('-hf', '--helix_file', dest = "helix_file", default=None, help='Input file for helices of your protein. To see the required format, check README or our GitHub page')
	args = parser.parse_args()

	###################################################################################################################
	assert len(args.output_name)>0, "Provide a name for output file"
	assert len(args.grofile)>0, "Provide input file"
	gro = MDAnalysis.Universe(args.grofile)
	list_of_non_ligands=["HOH","ARG","LYS","HIS","ASP","GLU","SER","THR", "ASN","GLN","PHE","TYR","TRP","CYS","GLY","PRO","ALA","VAL","ILE","LEU","MET"]
	potential_ligands={}
	i=0
	for resname in gro.residues.resnames:
	    if resname not in list_of_non_ligands and resname not in potential_ligands.values():
	        potential_ligands[i]=resname
	        i+=1
	for lig in potential_ligands:
	    print lig, potential_ligands[lig]
	ligand_name=potential_ligands[int(raw_input( "Choose a ligand to analyse:"))]


	if args.analysis_type=="occurance":
		md_sim = Topol_Data(args.grofile, None, ligand_name, args.offset)
		occurance = Occurance_analysis(args.grofile, args.xtcfile, ligand_name, args.cutoff, args.offset, md_sim)
	else: #if analysis type is anything different only one traj at time is going to be analysed
		assert args.xtcfile is None or len(args.xtcfile)<=1, "Only one trajectory at the time can be analysed."
		if args.xtcfile	is None:
			md_sim = Topol_Data(args.grofile, args.xtcfile, ligand_name, args.offset)
		else:
			md_sim = Topol_Data(args.grofile, args.xtcfile[0], ligand_name, args.offset)

		md_sim.find_res_to_plot(ligand_name, args.cutoff)
	md_sim.get_closest_ligand_atoms(ligand_name)

	hbonds = HBonds(md_sim, ligand_name)


	plots = Plots(md_sim)
	if args.diagram_type=="amino":
		plots.define_amino_acids()
		plots.plot_amino_diagramms()
	if args.diagram_type=="helices":
		assert len(args.helix_file)>0, "Provide a file defining helices"
		plots.define_helices(args.helix_file, args.offset)
		plots.plot_helices_diagramms()
	if args.diagram_type=="clock":
		assert args.analysis_type=="occurance", "Only occurance can be plotted with clock diagrams"
		plots.plot_clock_diagramms()


	molecule = Molecule(ligand_name, md_sim)
	if args.analysis_type=="RMSF" or args.analysis_type=="rmsf":
		rmsf = RMSF_measurements(md_sim)
		molecule = Molecule(ligand_name, md_sim, rmsf)
	molecule.convex_hull()
	molecule.make_new_projection_values()


	figure=Figure(molecule, hbonds)
	figure.add_bigger_box()
	figure.manage_the_plots()
	figure.draw_hbonds_in_graph()
	figure.put_everything_together()
	figure.write_final_draw_file(args.output_name)

	file_list=["molecule.svg"]
	for residue in md_sim.dict_of_plotted_res.keys():
		file_list.append(str(residue[3:])+".svg")
	for f in file_list:
		os.remove(f)

	print "Ready!"

	print "The figure you created has been saved under filename "+args.output_name+".svg"
	print "This file format can be converted to any other image file format using "
	print "free software such as ImageMagick, Inkscape or GIMP."
	print "    "
	print "#########################################################################"
	print "       "
	print "Thank you for using Ligand Interaction Network. Check our GitHub account"
	print "https://github.com/ldomic/lintools."
	
