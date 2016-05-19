import sys
reload(sys)
sys.setdefaultencoding('utf8')
import fileinput
import sys
import os
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import pylab
import numpy as np
from colors import *

class Figure(object):
    def __init__(self, molecule_object, diagram_type,topol_object, hbonds_object=None, plot_object=None, rmsf_object=None,tests=False):
        self.tests = tests
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
                    output1 = "<g transform='translate("+str(int(self.molecule.nearest_points_coords[residue][0]+(self.molecule.x_dim-900)/2)+46)+","+str(int(self.molecule.nearest_points_coords[residue][1]+(self.molecule.y_dim-450)/2)-54)+")'>" # add 100 to the left to have better alignment
                else:
                    output1 = "<g transform='translate("+str(int(self.molecule.nearest_points_coords[residue][0]+(self.molecule.x_dim-900)/2)-54)+","+str(int(self.molecule.nearest_points_coords[residue][1]+(self.molecule.y_dim-450)/2)-54)+")'>"
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
                    transform = "<g transform='translate("+str(int(self.molecule.nearest_points_coords[residue][0]+(self.molecule.x_dim-900)/2)-54)+","+str(int(self.molecule.nearest_points_coords[residue][1]+(self.molecule.y_dim-450)/2)-54)+")'>"  
                else:
                    transform = "<g transform='translate("+str(int(self.molecule.nearest_points_coords[residue][0]+(self.molecule.x_dim-900)/2)-54)+","+str(int(self.molecule.nearest_points_coords[residue][1]+(self.molecule.y_dim-450)/2)-54)+")'>"
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
            self.x_dim=self.molecule.x_dim
            self.molecule.x_dim=self.molecule.x_dim+100
        if diagram_type=="domains":
            self.x_dim=self.molecule.x_dim
            self.molecule.x_dim=self.molecule.x_dim+300
        if diagram_type=="amino":
            self.x_dim=self.molecule.x_dim
            self.y_dim=self.molecule.y_dim+60
        if diagram_type=="clock":
            self.x_dim=self.molecule.x_dim
            self.y_dim=self.molecule.y_dim
        start1 = "width='900px' height='450px' >"
        start2 = "<rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='900' height='450' x='0' y='0'> </rect>"
        bigger_box ="width='"+str(self.molecule.x_dim)+"px' height='"+str(self.molecule.y_dim+60)+"px' > "
        if self.rmsf!=None:
            big_box2= "<rect style='opacity:1.0;fill:white;stroke:none' width='"+str(self.molecule.x_dim)+"px' height='"+str(self.molecule.y_dim+60)+"px' x='0' y='0'> </rect> <g transform='translate("+str((self.x_dim-900)/2+100)+","+str((self.molecule.y_dim-450)/2)+")'>'<rect style='opacity:1.0;fill:#ffffff;stroke:none' width='900' height='450' x='0' y='0' /> "
        else:
            big_box2= "<rect style='opacity:1.0;fill:white;stroke:none' width='"+str(self.molecule.x_dim)+"px' height='"+str(self.molecule.y_dim+60)+"px' x='0' y='0'> </rect> <g transform='translate("+str((self.x_dim-900)/2)+","+str((self.molecule.y_dim-450)/2)+")'>'<rect style='opacity:1.0;fill:#ffffff;stroke:none' width='900' height='450' x='0' y='0' /> "
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
            if self.tests== False:
                coord=sys.argv[0][0:-11]+"legends/amino_legend.svg"
            else:
                coord=os.getcwd()+"/lintools/legends/amino_legend.svg"
            with open(coord,"r") as f:
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
            print float("{0:.1f}".format(self.rmsf.min_value))
            norm = matplotlib.colors.Normalize(vmin=float("{0:.1f}".format(self.rmsf.min_value)), vmax=float("{0:.1f}".format(self.rmsf.max_value)))
            cb1 = matplotlib.colorbar.ColorbarBase(ax1, cmap=new_cmap,
                                            norm=norm,
                                            orientation='vertical')
            cb1.ax.tick_params(labelsize=20) 
            cb1.set_label('RMSF',size=24, fontweight='bold')
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
         for atom in self.molecule.nearest_points_coords:
            self.white_circles = self.white_circles+"<circle cx='"+str(int(self.molecule.nearest_points_coords[atom][0])+(self.molecule.x_dim-900)/2)+"' cy='"+str(int(self.molecule.nearest_points_coords[atom][1])+(self.molecule.y_dim-450)/2)+"' r='30' fill='white' />"
    def draw_lines_in_graph(self):
        for residue in self.molecule.nearest_points_coords:
            self.draw_lines=self.draw_lines+"<line x1='"+str(int(self.molecule.nearest_points_coords[residue][0]))+"' y1='"+str(int(self.molecule.nearest_points_coords[residue][1]))+"' x2='"+str(float(self.molecule.ligand_atom_coords_from_diagr[self.topol.closest_atoms[residue][0][0]][0]))+"' y2='"+str(float(self.molecule.ligand_atom_coords_from_diagr[self.topol.closest_atoms[residue][0][0]][1]))+"' style='stroke:red;stroke-width:2' />"
    def draw_hbonds_in_graph(self):
        for bond in self.hbonds.hbonds_for_drawing:
            if bond[2]=="backbone":
                self.draw_lines=self.draw_lines+"<line x1='"+str(int(self.molecule.nearest_points_coords[bond[1]][0]))+"' y1='"+str(int(self.molecule.nearest_points_coords[bond[1]][1]))+"' x2='"+str(float(self.molecule.ligand_atom_coords_from_diagr[bond[0]][0]))+"' y2='"+str(float(self.molecule.ligand_atom_coords_from_diagr[bond[0]][1]))+"' style='stroke:blue;stroke-width:4' />"
            else:
                self.draw_lines=self.draw_lines+"<line stroke-dasharray='5,5'  x1='"+str(int(self.molecule.nearest_points_coords[bond[1]][0]))+"' y1='"+str(int(self.molecule.nearest_points_coords[bond[1]][1]))+"' x2='"+str(float(self.molecule.ligand_atom_coords_from_diagr[bond[0]][0]))+"' y2='"+str(float(self.molecule.ligand_atom_coords_from_diagr[bond[0]][1]))+"' style='stroke:blue;stroke-width:4' />"

    def put_everything_together(self):
        molecule_list = [self.filestart]+[self.draw_lines]+[self.draw_molecule]+[self.white_circles]+[self.draw_plots]+[self.legend]+[self.end_symbol]
        self.final_molecule = "".join(map(str,molecule_list))
    def write_final_draw_file(self, output_name):
        finalsvg = open(output_name+".svg","w")
        finalsvg.writelines(self.final_molecule)
        finalsvg.close()


class Residue_Info(object):
    def __init__(self, topol_object,occurrence, figure):
        self.residue_info = None
        self.universe = topol_object
        self.occurrence = occurrence
        self.figure = figure
        self.legend = ""
        self.gather_information()
        self.clean_up()
    def gather_information(self):
        self.plot_residue_on_off_rates()
        for residue in self.universe.dict_of_plotted_res:
            self.start_drawing(residue)
            self.write_final_draw_file(residue)
    def plot_residue_on_off_rates(self):
        self.colors = EarthSea[len(self.occurrence.residue_counts_on_off)]
        for res in self.universe.dict_of_plotted_res:
            for traj in self.occurrence.residue_counts_on_off:
                frame_count = int(self.occurrence.residue_counts_on_off[traj]["frame_count"])
                try:
                    Z=np.array([self.occurrence.residue_counts_on_off[traj][res]])
                except KeyError:
                    continue
                G = np.zeros((1,frame_count,3))
                G[Z>0.6] = [float(x)/256 for x in self.colors[traj-1][0]]
                G[Z==0.5] = [float(x)/256 for x in self.colors[traj-1][1]]
                G[Z<0.5] = [float(x)/256 for x in self.colors[traj-1][2]]
                fig, ax = plt.subplots(figsize=(4,0.35))
                ax = plt.imshow(G,interpolation='nearest',aspect="auto")
                ax.axes.get_yaxis().set_ticks([])
                if traj==1:
                    ax.axes.set_title(res, fontsize = 20)
                percent =str(int(float(self.occurrence.residue_counts[traj][res])/float(self.occurrence.residue_counts_on_off[traj]["frame_count"])*100))+"%"
                ax.axes.yaxis.set_label_position("right")
                ax.axes.set_ylabel(percent, rotation = "horizontal",fontsize = 20, labelpad = 35, va= "center",ha="center")
                div = int(len(ax.axes.get_xaxis().get_majorticklabels()))-3
                ns = int(self.occurrence.residue_counts_on_off[traj]["lastframe_time"])/1000
                step_by = ns/div
                end_value = ns+step_by
                ns_labels=[x for x in range(0,end_value,step_by)] 
                frame_count = int(self.occurrence.residue_counts_on_off[traj]["frame_count"])
                step_by_exist = frame_count/div
                existing_labels = [x for x in range(0,frame_count, step_by_exist)]
                plt.xticks(existing_labels,ns_labels)
                plt.axis("off")
                pylab.savefig(str(res)+"_"+str(traj)+".svg")

    def start_drawing(self, residue):
        file_width = 320
        self.end_symbol = "</svg>"
        self.on_off_plots = ""
        resname_svg = "<text style='font-family:Times New Roman;font-size:24.0px;' x='45' y='20'>" +str(residue)+  "</text>"
        
        y=25
        y_for_text = 42
        for traj in self.occurrence.residue_counts_on_off:
            try:
                with open(str(residue)+"_"+str(traj)+".svg", "r") as f:
                    lines = f.readlines()
                    self.on_off_plots = self.on_off_plots+ "<g  transform='translate(0,"+str(y)+")'>\n"+"".join(map(str,lines[21:23]))+"\n</g>\n"
                    percent =str(int(float(self.occurrence.residue_counts[traj][residue])/float(self.occurrence.residue_counts_on_off[traj]["frame_count"])*100))+"%"
                    self.on_off_plots = self.on_off_plots + "<text style='font-family:Monospace;font-size:14.0px;' x='260' y='"+str(y_for_text)+"'> "+percent+" </text>"
                    f.close()
                y+=25
                y_for_text+=25
            except IOError:
                continue
        #make a legend
        self.legend = "<g transform='translate(0,"+str(y+10)+")'> "
        x=35
        for traj in self.occurrence.residue_counts_on_off:
            self.legend = self.legend+"<rect style='opacity:1.0;fill:rgb("+",".join(str(i) for i in self.colors[traj-1][0])+");stroke:none' width='15' height='15' x='"+str(x)+"' y='0' /> "
            x+=15
        self.legend = self.legend + "<text style='font-family:Monospace;font-size:10.0px;' x='"+str(x)+"' y='10'>"+"&lt;"+str(self.occurrence.residue_counts_on_off[traj]["cutoff"])+u"\u212B"+"</text>"
        self.legend = self.legend + "</g>"
        x+=40
        self.legend = self.legend +"<g transform='translate(0,"+str(y+10)+")'> "
        for traj in self.occurrence.residue_counts_on_off:
            self.legend = self.legend+"<rect style='opacity:1.0;fill:rgb("+",".join(str(i) for i in self.colors[traj-1][1])+");stroke:none' width='15' height='15' x='"+str(x)+"' y='0' /> "
            x+=15
        self.legend = self.legend + "<text style='font-family:Monospace;font-size:10.0px;' x='"+str(x-5)+"' y='10'> "+str(self.occurrence.residue_counts_on_off[traj]["cutoff"])+"-"+str(self.occurrence.residue_counts_on_off[traj]["cutoff"]+self.occurrence.residue_counts_on_off[traj]["cutoff_expanded"])+u"\u212B"+"</text>"
        self.legend = self.legend + "</g>"

        self.filestart = "<svg version='1.1' baseProfile='full' xmlns:svg='http://www.w3.org/2000/svg' xmlns:rdkit='http://www.rdkit.org/xml' xmlns:xlink='http://www.w3.org/1999/xlink' xml:space='preserve' width='"+str(file_width)+"' height='"+str(y+60)+"' > "
        big_box = "<rect style='opacity:1.0;fill:#ffffff;stroke:black;stroke-width:3' width='280' height='"+str(y+20)+"' x='25' y='10' /> "
        small_box = "<rect style='opacity:1.0;fill:#ffffff;stroke:white;stroke-width:3' width='95' height='30' x='40' y='0' /> "
        self.filestart = self.filestart + big_box + small_box + resname_svg
        #deal with 2D molecule
        #self.molecule = "<g  transform='scale(0.3) translate(0,"+str(y+225)+")'>\n"+ self.figure.draw_molecule

    def write_final_draw_file(self, residue):
        file_list = [self.filestart]+[self.on_off_plots]+[self.legend]+[self.end_symbol]
        self.final_residue_file = "".join(map(str,file_list))
        finalsvg = open(str(residue)+".svg","w")
        finalsvg.writelines(self.final_residue_file)
        finalsvg.close()
        

    def clean_up(self):
        for res in self.universe.dict_of_plotted_res:
            for traj in self.occurrence.residue_counts_on_off:
                if os.path.isfile(str(res)+"_"+str(traj)+".svg")==True:
                    os.remove(str(res)+"_"+str(traj)+".svg")

