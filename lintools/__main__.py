import sys
import lintools
import MDAnalysis
from argparse import ArgumentParser


def main(args=None):
    """The main routine."""
    if args is None:
        args = sys.argv[1:]

    print("This is the main routine.")
    print("It should do something interesting.")
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
if __name__ == "__main__":
    main()
