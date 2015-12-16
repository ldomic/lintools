
import fileinput
import sys

class Figure(object):
    def __init__(self, molecule_object, diagram_type, hbonds_object=None, plot_object=None):
        self.draw_plots = None
        self.draw_molecule =None
        self.draw_lines=" "
        self.final_molecule =None
        self.molecule = molecule_object
        self.hbonds = hbonds_object
        self.plots = plot_object
        self.legend = None
        self.make_legends(diagram_type)
        self.add_bigger_box(diagram_type)
    def change_lines_in_svg(self,filename, string1,string2):
        for i,line in enumerate(fileinput.input(filename, inplace=1)):
            sys.stdout.write(line.replace(str(string1),str(string2)))
    def manage_the_plots(self):
        diagram=""
        for residue in self.molecule.nearest_points_coords.keys():
            for i, line in enumerate(fileinput.input(str(residue[3:])+".svg", inplace=1)):
                if i <= 8:
                    continue
                else:
                    sys.stdout.write(line.replace ("</svg>","</g>"))
            input1 = "</defs>"
            output1 = "<g transform='translate("+str(int(self.molecule.nearest_points_coords[residue][0]+(self.molecule.x_dim-600)/2)-54)+","+str(int(self.molecule.nearest_points_coords[residue][1]+(self.molecule.y_dim-300)/2)-54)+")'>"
            self.change_lines_in_svg(str(residue[3:])+'.svg', input1, output1)
            input2 = "font-style:normal;"
            output2 = "font-style:normal;font-weight:bold;"
            self.change_lines_in_svg(str(residue[3:])+'.svg', input2, output2)
            with open(str(residue[3:])+".svg", "r") as f:
                lines = f.readlines()
                diagram = diagram +"".join(map(str,lines))
                f.close()
        self.draw_plots = diagram
    def add_bigger_box(self, diagram_type):
        """Rewrite the molecule.svg file line by line, otherwise it fails."""
        start1 = "width='600px' height='300px' >"
        start2 = "<rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='600' height='300' x='0' y='0'> </rect>"
        bigger_box ="width='"+str(self.molecule.x_dim)+"px' height='"+str(self.molecule.y_dim)+"px' > "
        big_box2= "<rect style='opacity:1.0;fill:white;stroke:none' width='"+str(self.molecule.x_dim)+"px' height='"+str(self.molecule.y_dim)+"px' x='0' y='0'> </rect><g transform='translate("+str((self.molecule.x_dim-600)/2)+","+str((self.molecule.y_dim-300)/2)+")'>'<rect style='opacity:1.0;fill:#ffffff;stroke:none' width='600' height='300' x='0' y='0'> </rect>"
        domain_box = "width='"+str(self.molecule.x_dim+300)+"px' height='"+str(self.molecule.y_dim)+"px' > "
        amino_box = "width='"+str(self.molecule.x_dim)+"px' height='"+str(self.molecule.y_dim+60)+"px' > "
        self.end_symbol = "</svg>"
        no_end_symbol = "</g>"
        num_lines = sum(1 for line in open('molecule.svg'))
        if diagram_type=="domains":
            self.change_lines_in_svg("molecule.svg", start1, domain_box)
        if diagram_type=="amino":
            self.change_lines_in_svg("molecule.svg", start1, amino_box)
        else:
            self.change_lines_in_svg("molecule.svg", start1, bigger_box)
        self.change_lines_in_svg("molecule.svg", start2, big_box2)
        self.change_lines_in_svg("molecule.svg", self.end_symbol, no_end_symbol)
        with open("molecule.svg","r") as f:
            lines = f.readlines()
            self.draw_molecule ="".join(map(str,lines))
            f.close()
    def make_legends(self, diagram_type,domain_file=None):
        if diagram_type=="amino":
            self.legend = "<g transform='translate(0,"+str(self.molecule.y_dim+20)+")'>"
            with open("legends/amino_legend.svg","r") as f:
                lines = f.readlines()
                legend ="".join(map(str,lines))
                f.close()
            self.legend = self.legend+legend
        if diagram_type=="domains":
            print "In progress"
            sorted_dom = sorted(self.plots.plotted_domains)
            #continue here... The code does not work - try HSG for colors
            #self.legend="<g>"
            self.legend="<g transform='translate("+str(self.molecule.x_dim)+",0)'>"
            y=50
            for dom in sorted_dom:
                self.legend=self.legend+"<circle cx='25' cy='"+str(y)+"' r='20' stroke='"+str(dom[2])+"' stroke-width='5' fill='none' />"+"  <text x='50' y='"+str(y+5)+"' style='font-size:20px;fill:#000000;font-family:Verdana'>"+str(dom[1])+"</text>"
                y+=50
            self.legend=self.legend+"</g>"
    def draw_lines_in_graph(self):
        for residue in self.molecule.nearest_points_coords:
            self.draw_lines=self.draw_lines+"<line x1='"+str(int(self.molecule.nearest_points_coords[residue][0])+(self.molecule.x_dim-600)/2)+"' y1='"+str(int(self.molecule.nearest_points_coords[residue][1])+(self.molecule.y_dim-300)/2)+"' x2='"+str(float(self.molecule.atom_coords_from_diagramm[residue][0])+(self.molecule.x_dim-600)/2)+"' y2='"+str(float(self.molecule.atom_coords_from_diagramm[residue][1])+(self.molecule.y_dim-300)/2)+"' style='stroke:'red';stroke-width:2' />"
    def draw_hbonds_in_graph(self):
        for bond in self.hbonds.hbonds_for_drawing:
            self.draw_lines=self.draw_lines+"<line stroke-dasharray='5,5'  x1='"+str(int(self.molecule.nearest_points_coords[bond[1]][0])+(self.molecule.x_dim-600)/2)+"' y1='"+str(int(self.molecule.nearest_points_coords[bond[1]][1])+(self.molecule.y_dim-300)/2)+"' x2='"+str(float(self.molecule.atom_coords_from_diagramm[bond[1]][0])+(self.molecule.x_dim-600)/2)+"' y2='"+str(float(self.molecule.atom_coords_from_diagramm[bond[1]][1])+(self.molecule.y_dim-300)/2)+"' style='stroke:black;stroke-width:4' />"
    def put_everything_together(self):
        molecule_list = [self.draw_molecule]+[self.draw_lines]+[self.draw_plots]+[self.legend]+[self.end_symbol]
        self.final_molecule = "".join(map(str,molecule_list))
    def write_final_draw_file(self, output_name):
        finalsvg = open(output_name+".svg","w")
        finalsvg.writelines(self.final_molecule)
        finalsvg.close()