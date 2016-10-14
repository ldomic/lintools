import fileinput
import sys
class Figure(object):
    """
    This class inherits Molecule(),HBonds() and Data() class objects that provide the information necessary
    to put together the final SVG image. This is done by compiling and altering the existing SVG information
    as well as creating new information (e.g. lines for hydrogen bonds).
    """
    __version__="09.2016"
    def __init__(self,molecule_object,topology_data_object,hydrogen_bonds_object):
        self.molecule = molecule_object
        self.topology_data = topology_data_object
        self.hbonds = hydrogen_bonds_object
        self.draw_plots = None
        self.draw_lines = ""
        self.white_circles = ""
    def change_lines_in_svg(self,filename, string1,string2):
        """Used to change lines in an SVG file. String1 is input - what is already present
        in the file, while string2 is output - the desired line. The input has to be searchable
        since .replace function is used. 
        Takes:
            * filename * - name of the SVG file to change
            * string1 * - input string
            * string2 * - output string

        """
        for i,line in enumerate(fileinput.input(filename, inplace=1)):
            sys.stdout.write(line.replace(str(string1),str(string2)))
    def manage_the_plots(self):
        """
        Each plotted residue SVG file is edited and the contents of the file transfered into 
        an SVG group (g). This allows a plot to be included in the final image and the group
        transform function allows to move the plot to the correct 2D coordinates in respect 
        to the 2D drug molecule. The groups are added to self.draw_plots variable which is 
        included in the final figure.
        """
        diagram = ""
        for residue in sorted(self.molecule.nearest_points_coords.keys()):
            for i, line in enumerate(fileinput.input(residue[1]+residue[2]+".svg", inplace=1)):
                    if i <= 8:
                        continue
                    else:
                        sys.stdout.write(line.replace ("</svg>","</g>"))
            input1 = "</defs>"
            output1 = "<g transform='translate("+str(int(self.molecule.nearest_points_coords[residue][0]+(self.molecule.x_dim-self.molecule.molsize1)/2)-90)+","+str(int(self.molecule.nearest_points_coords[residue][1]+(self.molecule.y_dim-self.molecule.molsize2)/2)-90)+")'>"
            self.change_lines_in_svg(residue[1]+residue[2]+'.svg', input1, output1)
            input2 = "font-style:normal;"
            output2 = "font-style:normal;font-weight:bold;"
            self.change_lines_in_svg(residue[1]+residue[2]+'.svg', input2, output2)
            with open(residue[1]+residue[2]+".svg", "r") as f:
                lines = f.readlines()
                diagram = diagram +"".join(map(str,lines))
                f.close()
        self.draw_plots = diagram
    def add_bigger_box(self):
        """
        Sets the size of the figure by expanding the space of molecule.svg file. These dimension have been 
        previously determined. Also makes the lines of the molecule thicker.
        """
        start1 = "width='"+str(self.molecule.molsize1)+"px' height='"+str(self.molecule.molsize2)+"px' >"
        start2 = "<rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='"+str(int(self.molecule.molsize1))+"' height='"+str(int(self.molecule.molsize2))+"' x='0' y='0'> </rect>"
        bigger_box ="width='"+str(int(self.molecule.x_dim))+"px' height='"+str(int(self.molecule.y_dim))+"px' > "
        big_box2= "<rect style='opacity:1.0;fill:white;stroke:none' width='"+str(int(self.molecule.x_dim))+"px' height='"+str(int(self.molecule.y_dim))+"px' x='0' y='0'> </rect> <g transform='translate("+str((self.molecule.x_dim-self.molecule.molsize1)/2)+","+str((self.molecule.y_dim-self.molecule.molsize2)/2)+")'>'<rect style='opacity:1.0;fill:#ffffff;stroke:none' width='"+str(self.molecule.molsize1)+"' height='"+str(self.molecule.molsize2)+"' x='0' y='0' /> "
        self.end_symbol = "</svg>"
        no_end_symbol = "</g>"
        #Make the lines in molecule drawing thicker to look better with the large plots
        linewidth1 = "stroke-width:2px"
        linewidth2 = "stroke-width:5px"
        self.change_lines_in_svg("molecule.svg", linewidth1,linewidth2)

        self.change_lines_in_svg("molecule.svg", start1, bigger_box)
        self.change_lines_in_svg("molecule.svg", start2, big_box2)
        self.change_lines_in_svg("molecule.svg", self.end_symbol, no_end_symbol)
        with open("molecule.svg","r") as f:
            lines = f.readlines()
            self.filestart = " ".join(map(str,lines[0:8]))
            self.draw_molecule ="".join(map(str,lines[8:]))
            f.close()
    def draw_white_circles(self):
        """
        The plots are transparent and therefore this function is required to cover up the middle 
        part of the plots (that contains text). This function creates white circles that provide
        background.
        """
        for atom in sorted(self.molecule.nearest_points_coords.keys()):
                self.white_circles = self.white_circles+"<circle cx='"+str(int(self.molecule.nearest_points_coords[atom][0]))+"' cy='"+str(int(self.molecule.nearest_points_coords[atom][1]))+"' r='55' fill='white' />"
    def draw_hydrogen_bonds(self):
        """
        For each bond that has been determined to be important, a line gets drawn.
        """
        if self.hbonds!=None:
            for bond in self.hbonds.hbonds_for_drawing:
                atom = self.topology_data.universe.atoms[bond[0]-1] #zero-based index vs one-based index
                residue = (atom.resname, str(atom.resid), atom.segid)
                if bond[2] in ["N","O","H"]:
                    #backbone interactions
                    self.draw_lines=self.draw_lines+"<line x1='"+str(int(self.molecule.nearest_points_coords[residue][0]))+"' y1='"+str(int(self.molecule.nearest_points_coords[residue][1]))+"' x2='"+str(float(self.molecule.ligand_atom_coords_from_diagr[bond[1]][0]))+"' y2='"+str(float(self.molecule.ligand_atom_coords_from_diagr[bond[1]][1]))+"' style='stroke:black;stroke-width:4' />"
                else:    
                    #sidechain interactions
                    self.draw_lines=self.draw_lines+"<line stroke-dasharray='5,5'  x1='"+str(int(self.molecule.nearest_points_coords[residue][0]))+"' y1='"+str(int(self.molecule.nearest_points_coords[residue][1]))+"' x2='"+str(float(self.molecule.ligand_atom_coords_from_diagr[bond[1]][0]))+"' y2='"+str(float(self.molecule.ligand_atom_coords_from_diagr[bond[1]][1]))+"' style='stroke:black;stroke-width:4' />"

    def put_everything_together(self):
        """
        All of the elements of the final SVG file are put together in the correct order (e.g. lines are placed behind plots 
        and the molecule).
        """
        molecule_list = [self.filestart]+[self.draw_lines]+[self.white_circles]+[self.draw_molecule]+[self.draw_plots]+[self.end_symbol]
        self.final_molecule = "".join(map(str,molecule_list))
    def write_final_draw_file(self, output_name):
        """The result of put_everything_together() function is writen in a file.
        Takes:
            * output_name * - name for the output file 
        """
        finalsvg = open(output_name+".svg","w")
        finalsvg.writelines(self.final_molecule)
        finalsvg.close()