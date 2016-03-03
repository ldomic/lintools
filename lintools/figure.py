
import fileinput
import sys
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import pylab
import numpy as np

class Figure(object):
    def __init__(self, molecule_object, diagram_type,topol_object, hbonds_object=None, plot_object=None, rmsf_object=None):
        self.draw_plots = None
        self.draw_molecule =None
        self.draw_lines=" "
        self.filestart=""
        self.white_circles=""
        self.final_molecule =None
        self.topol = topol_object
        self.molecule = molecule_object
        self.hbonds = hbonds_object
        self.plots = plot_object
        self.rmsf = rmsf_object
        self.legend = ""
        self.make_legends(diagram_type)
        self.manage_the_plots(diagram_type)
        self.add_bigger_box(diagram_type)
    def change_lines_in_svg(self,filename, string1,string2):
        for i,line in enumerate(fileinput.input(filename, inplace=1)):
            sys.stdout.write(line.replace(str(string1),str(string2)))
    def manage_the_plots(self, diagram_type):
        diagram=""
        if diagram_type!="domains":
            for residue in self.molecule.nearest_points_coords.keys():
                for i, line in enumerate(fileinput.input(str(residue[3:])+".svg", inplace=1)):
                    if i <= 8:
                        continue
                    else:
                        sys.stdout.write(line.replace ("</svg>","</g>"))
                input1 = "</defs>"
                if self.rmsf!=None:
                    output1 = "<g transform='translate("+str(int(self.molecule.nearest_points_coords[residue][0]+(self.molecule.x_dim-600)/2)+46)+","+str(int(self.molecule.nearest_points_coords[residue][1]+(self.molecule.y_dim-300)/2)-54)+")'>" # add 100 to the left to have better alignment
                else:
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
        else:
            self.draw_plots=""
            for residue in self.plots.residues_within_domain:
                if self.rmsf!=None:
                    transform = "<g transform='translate("+str(int(self.molecule.nearest_points_coords[residue][0]+(self.molecule.x_dim-600)/2)+46)+","+str(int(self.molecule.nearest_points_coords[residue][1]+(self.molecule.y_dim-300)/2)-54)+")'>"  
                else:
                    transform = "<g transform='translate("+str(int(self.molecule.nearest_points_coords[residue][0]+(self.molecule.x_dim-600)/2)-54)+","+str(int(self.molecule.nearest_points_coords[residue][1]+(self.molecule.y_dim-300)/2)-54)+")'>"
                path = "<rect style='fill:none' width='108' height='108' x='0' y='0' />"
                circle = "<circle cx='54' cy='54' r='38' stroke='"+str(self.plots.residues_within_domain[residue][2])+"' stroke-width='10' fill='white' />"
                dashed_circle =  "<circle cx='54' cy='54' r='38' stroke='"+str(self.plots.residues_within_domain[residue][2])+"' stroke-width='10' fill='white' stroke-dasharray='10,5'" 
                if str(self.plots.residues_within_domain[residue][3][0][0])=="Y":
                    self.draw_plots = self.draw_plots+transform+path+dashed_circle
                else:
                    self.draw_plots = self.draw_plots+transform+path+circle
                with open(str(residue[3:])+".svg", "r") as f:
                    lines=f.readlines()
                    for line in lines:
                        if line.startswith("    <text"):
                            self.draw_plots=self.draw_plots+line
                self.draw_plots=self.draw_plots+"</g>"
    def add_bigger_box(self, diagram_type):
        """Rewrite the molecule.svg file line by line, otherwise it fails."""
        if self.rmsf!=None:
            x_dim=self.molecule.x_dim
            self.molecule.x_dim=self.molecule.x_dim+100
        if diagram_type=="domains":
            x_dim=self.molecule.x_dim
            self.molecule.x_dim=self.molecule.x_dim+300
        if diagram_type=="amino":
            y_dim=self.molecule.y_dim+60
        start1 = "width='600px' height='300px' >"
        start2 = "<rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='600' height='300' x='0' y='0'> </rect>"
        bigger_box ="width='"+str(self.molecule.x_dim)+"px' height='"+str(self.molecule.y_dim)+"px' > "
        if self.rmsf!=None:
            big_box2= "<rect style='opacity:1.0;fill:white;stroke:none' width='"+str(self.molecule.x_dim)+"px' height='"+str(self.molecule.y_dim)+"px' x='0' y='0'> </rect> <g transform='translate("+str((x_dim-600)/2+100)+","+str((self.molecule.y_dim-300)/2)+")'>'<rect style='opacity:1.0;fill:#ffffff;stroke:none' width='600' height='300' x='0' y='0' /> "
        else:
            big_box2= "<rect style='opacity:1.0;fill:white;stroke:none' width='"+str(self.molecule.x_dim)+"px' height='"+str(self.molecule.y_dim)+"px' x='0' y='0'> </rect> <g transform='translate("+str((self.molecule.x_dim-600)/2)+","+str((self.molecule.y_dim-300)/2)+")'>'<rect style='opacity:1.0;fill:#ffffff;stroke:none' width='600' height='300' x='0' y='0' /> "
        self.end_symbol = "</svg>"
        no_end_symbol = "</g>"
        self.change_lines_in_svg("molecule.svg", start1, bigger_box)
        self.change_lines_in_svg("molecule.svg", start2, big_box2)
        self.change_lines_in_svg("molecule.svg", self.end_symbol, no_end_symbol)
        with open("molecule.svg","r") as f:
            lines = f.readlines()
            self.filestart = " ".join(map(str,lines[0:8]))
            self.draw_molecule ="".join(map(str,lines[8:]))
            f.close()
    def make_legends(self, diagram_type,domain_file=None):
        if diagram_type=="amino":
            self.legend = "<g transform='translate(0,"+str(self.molecule.y_dim+20)+")'>"
            with open(sys.argv[0][0:-11]+"legends/amino_legend.svg","r") as f:
                lines = f.readlines()
                legend ="".join(map(str,lines))
                f.close()
            self.legend = self.legend+legend
        if diagram_type=="domains":
            if self.rmsf!=None:
                x_dim=self.molecule.x_dim+100
                y_dim= (self.molecule.y_dim-519.12-56.88)/2
            else:
                x_dim=self.molecule.x_dim
                y_dim=0
            sorted_dom = sorted(self.plots.plotted_domains)
            self.legend="<g transform='translate("+str(x_dim)+","+str(y_dim)+")'>"
            y=50
            for dom in sorted_dom:
                if dom[3][0][0]=="Y":
                    self.legend=self.legend+"<circle cx='25' cy='"+str(y)+"' r='20' stroke='"+str(dom[2])+"' stroke-width='5' stroke-dasharray='10,5' fill='none' />"+"  <text x='50' y='"+str(y+5)+"' style='font-size:20px;fill:#000000;font-family:Verdana'>"+str(dom[1])+"</text>"
                else:
                    self.legend=self.legend+"<circle cx='25' cy='"+str(y)+"' r='20' stroke='"+str(dom[2])+"' stroke-width='5' fill='none' />"+"  <text x='50' y='"+str(y+5)+"' style='font-size:20px;fill:#000000;font-family:Verdana'>"+str(dom[1])+"</text>"
                y+=50
            self.legend=self.legend+"</g>"
        if self.rmsf!=None:
            print "Making colorbar"
            fig = plt.figure(figsize=(1, 8))
            ax1 = fig.add_axes([0, 0, 0.5, 0.9])
            cmap = plt.get_cmap("hsv_r")
            new_cmap = self.truncate_colormap(cmap,0.67,1.00)
            norm = matplotlib.colors.Normalize(vmin=int(self.rmsf.min_value), vmax=int(self.rmsf.max_value))
            cb1 = matplotlib.colorbar.ColorbarBase(ax1, cmap=new_cmap,
                                            norm=norm,
                                            orientation='vertical')
            cb1.set_label('RMSF',size=19, fontweight='bold')
            pylab.savefig("rmsf_colorbar.svg", dpi=100, transparent=True)
            self.legend=self.legend+self.manage_the_rmsf_colorbar()


    def truncate_colormap(self, cmap, minval=0.0, maxval=1.0, n=100):
        new_cmap = colors.LinearSegmentedColormap.from_list(
            'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
            cmap(np.linspace(minval, maxval, n)))
        return new_cmap
    def manage_the_rmsf_colorbar(self):
        for i, line in enumerate(fileinput.input("rmsf_colorbar.svg", inplace=1)):
            if i <= 19:
                continue
            else:
                sys.stdout.write(line.replace ("</svg>","</g>"))
        y_dim = (self.molecule.y_dim-519.12-56.88)/2
        with open("rmsf_colorbar.svg", "r") as f:
            lines = f.readlines()
            colorbar = "<g transform='translate(10,"+str(y_dim)+")'> "+"".join(map(str,lines))
            f.close()
        return colorbar
    def draw_white_circles_at_atoms(self):
        if self.rmsf!=None:
            for atom in self.molecule.nearest_points_coords:
                self.white_circles = self.white_circles+"<circle cx='"+str((int(self.molecule.nearest_points_coords[atom][0])+(self.molecule.x_dim-600)/2)+100)+"' cy='"+str(int(self.molecule.nearest_points_coords[atom][1])+(self.molecule.y_dim-300)/2)+"' r='30' fill='white' />"
        else:
             for atom in self.molecule.nearest_points_coords:
                self.white_circles = self.white_circles+"<circle cx='"+str(int(self.molecule.nearest_points_coords[atom][0])+(self.molecule.x_dim-600)/2)+"' cy='"+str(int(self.molecule.nearest_points_coords[atom][1])+(self.molecule.y_dim-300)/2)+"' r='30' fill='white' />"
    def draw_lines_in_graph(self):
        for residue in self.molecule.nearest_points_coords:
            self.draw_lines=self.draw_lines+"<line x1='"+str(int(self.molecule.nearest_points_coords[residue][0]))+"' y1='"+str(int(self.molecule.nearest_points_coords[residue][1]))+"' x2='"+str(float(self.molecule.atom_coords_from_diagramm[residue][0]))+"' y2='"+str(float(self.molecule.atom_coords_from_diagramm[residue][1]))+"' style='stroke:red;stroke-width:2' />"
    def draw_hbonds_in_graph(self):
        for bond in self.hbonds.hbonds_for_drawing:
            self.draw_lines=self.draw_lines+"<line stroke-dasharray='5,5'  x1='"+str(int(self.molecule.nearest_points_coords[bond[1]][0]))+"' y1='"+str(int(self.molecule.nearest_points_coords[bond[1]][1]))+"' x2='"+str(float(self.molecule.ligand_atom_coords_from_diagr[bond[0]][0]))+"' y2='"+str(float(self.molecule.ligand_atom_coords_from_diagr[bond[0]][1]))+"' style='stroke:black;stroke-width:4' />"
    def put_everything_together(self):
        molecule_list = [self.filestart]+[self.draw_lines]+[self.draw_molecule]+[self.white_circles]+[self.draw_plots]+[self.legend]+[self.end_symbol]
        self.final_molecule = "".join(map(str,molecule_list))
    def write_final_draw_file(self, output_name):
        finalsvg = open(output_name+".svg","w")
        finalsvg.writelines(self.final_molecule)
        finalsvg.close()