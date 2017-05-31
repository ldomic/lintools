import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import pylab
import numpy


class Plots(object):
    """
    This module plots the residue data.
    Types that are available are:

    * amino * (09.2016) - circle with residue name and id (and chain) colored by residue type
    * domain * (09.2016) - residue name and id (and chain) circled by differently colored rings,
    ring color depends on the chain the residue is on.
    * clock * (09.2016) - residue name and id (and chain) circled by one or several differently
    colored rings that display the fraction of time this residue has spent in the vicinity of
    the ligand.

    Takes:
        * topology_data_object * - information about the system (lintools.Data object)
        * diagram_type * - string of the selected diagram type (e.g. "amino" or "clocks")
    """
    __version__ = "09.2016"
    matplotlib.rcParams['svg.fonttype'] = 'none'
    matplotlib.rcParams['font.weight']=900
    matplotlib.rcParams['text.usetex'] = False
    matplotlib.rcParams['patch.linewidth'] = 0
    def __init__(self, topology_data_object,diagram_type,colormap='summer'):
        self.topology_data = topology_data_object
        self.colors_amino_acids = {"acidic":"#D9774B", "basic":"#889DCC",
                                   "aromatic":"#9FC74A", "polar":"#D06AC1",
                                   "hydrophobic":"#6AC297","lipids":"#ffff99",
                                   "water":"turquoise","ions":"gold"}
        self.amino_acids = {"ASP":"acidic","GLU":"acidic","LYS":"basic","ARG":"basic",
                       "PHE":"aromatic","TYR":"aromatic","TRP":"aromatic","SER":"polar",
                       "THR":"polar","ASN":"polar","GLN":"polar","CYS":"polar",
                       "HIS":"polar","ALA":"hydrophobic","VAL":"hydrophobic",
                       "ILE":"hydrophobic","LEU":"hydrophobic","MET":"hydrophobic","GLY":"hydrophobic","PRO":"hydrophobic",
                       "PC":"lipids","HOH":"water","SOL":"water"}
        if diagram_type == "amino":
            self.plot_amino_diagrams()
        if diagram_type == "domains":
            self.plot_domain_diagrams()
        if diagram_type == "clock":
            self.plot_clock_diagrams(colormap)
    def plot_amino_diagrams(self):
        """
        Plotting of amino diagrams - circles with residue name and id, colored according to the
        residue type. If the protein has more than one chain, chain identity is also included in
        the plot. The plot is saved as svg file with residue id and chain id as filename for more
        certain identification.
        """

        for res in self.topology_data.dict_of_plotted_res:
            try:
                color = [self.colors_amino_acids[self.amino_acids[res[0]]],'white']
            except KeyError:
                color = ["pink",'white']
            plt.figure(figsize=(2.5,2.5))
            ring1,_=plt.pie([1],  radius=1, startangle=90, colors=color, counterclock=False)
            plt.axis('equal')
            plt.setp(ring1, width=1, edgecolor=color[0])
            if len(self.topology_data.universe.protein.segments)<=1:
                #Parameters for amino diagrams without segids
                plt.text(0,-0.45,res[0]+"\n"+res[1],ha='center',size=36, fontweight="bold")
            else:
                #Parameters for amino diagrams with segids
                plt.text(0,-0.37,res[0]+"\n"+res[1]+" "+res[2],ha='center',size=30, fontweight="bold")
            #play with the dpi
            pylab.savefig(str(res[1])+res[2]+".svg", dpi=300, transparent=True)


    def plot_domain_diagrams(self):
        """
        Plotting domain diagrams - a ring around the residue name and id and chain id. The colors are
        determined by the chain id and are extracted from matplotlib colormap "terrain" (ver. "09.2016")
        The plot is saved as svg file with residue id and chain id as filename for more certain
        identification.
        """
        # width of the circle around plot
        width=0.20


        # define color library
        cmap = plt.get_cmap('terrain')
        colors = [cmap(i) for i in numpy.linspace(0, 0.75, len(self.topology_data.universe.protein.segments))]
        domain_colors = {seg:colors[i] for i,seg in enumerate(self.topology_data.universe.protein.segments.segids.tolist())}
        for res in self.topology_data.dict_of_plotted_res:
            color = [domain_colors[res[2]],'white']
            #color = [self.colors_amino_acids[self.amino_acids[res.resname]],'white']
            plt.figure(figsize=(2.5,2.5))
            ring1,_=plt.pie([1],  radius=1-width, startangle=90, colors=color, counterclock=False)
            plt.axis('equal')
            plt.setp(ring1, width=width, edgecolor='white')
            if len(self.topology_data.universe.protein.segments)<=1:
                #Parameters for amino diagrams without segids
                plt.text(0,-0.4,res[0]+"\n"+res[1],ha='center',size=36, fontweight="bold")
            else:
                #Parameters for amino diagrams with segids
                plt.text(0,-0.22,res[0]+"\n"+res[1]+"\n"+res[2],ha='center',size=28, fontweight="bold")
            #play with the dpi
            pylab.savefig(res[1]+res[2]+".svg", dpi=300, transparent=True)


    def plot_clock_diagrams(self,colormap):
        """
        Ploting clock diagrams - one or more rings around residue name and id (and chain id).
        The rings show the fraction of simulation time this residue has spent in the vicinity of the
        ligand - characterised by distance.
        """
        cmap = plt.get_cmap(colormap)
        for res in self.topology_data.dict_of_plotted_res:
            colors = [cmap(i) for i in numpy.linspace(0, 1, len(self.topology_data.dict_of_plotted_res[res]))]
            traj_colors_ = {traj:colors[i] for i,traj in enumerate(self.topology_data.dict_of_plotted_res[res])}
            plt.figure(figsize=(2.25, 2.25))
            ring_number=[sum(1 for x in v if x) for k,v in self.topology_data.dict_of_plotted_res.items()][0]
            self.topology_data.ring_number = ring_number
            rings=[]
            # When only a few rings to plot they can be thicker
            if ring_number<2:
                width = 0.3
            else:
                width = 0.2
            for ring in range(0,ring_number):
                ring,_=plt.pie([self.topology_data.dict_of_plotted_res[res][ring],1-self.topology_data.dict_of_plotted_res[res][ring]],  radius=0.9+width*(ring+1), startangle=90, colors=[colors[ring],"white"], counterclock=False)
                rings=rings+ring
            plt.setp(rings, width=width)
            if len(self.topology_data.universe.protein.segments)<=1:
            #Settings with domain
                plt.text(-0.0,-0.62,res[0]+"\n"+res[1],ha='center',size=32, fontweight='bold')
            else:
                plt.text(-0.0,-0.72,res[0]+"\n"+res[1]+"\n"+res[2],ha='center',size=25, fontweight='bold')
            pylab.savefig(res[1]+res[2]+".svg", dpi=300, transparent=True)
