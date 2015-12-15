#sudo apt-get install -y python-qt4
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import pylab
from topol import Topol_Data
class Plots(object):
    matplotlib.rcParams['svg.fonttype'] = 'none'
    matplotlib.rcParams['font.weight']=900
    matplotlib.rcParams['text.usetex'] = False
    matplotlib.rcParams['patch.linewidth'] = 0  
    def __init__(self, topol_object):
        #Group amino acids by charge and structure
        self.amino_acids = {"acidic":["ASP","GLU"], "basic":["LYS","ARG"], "aromatic":["PHE","TYR","TRP"],"polar":["SER","THR","ASN","GLN","CYS","HIS"],"hydrophobic":["ALA","VAL","ILE","LEU","MET","GLY"]}
        self.colors_amino_acids = {"acidic":"#D9774B", "basic":"#889DCC", "aromatic":"#9FC74A", "polar":"#D06AC1","hydrophobic":"#6AC297"}
        self.amino_acid_type={}
        self.colors_helices ={1:(0.5490,0.7961, 0.8),2:(0.7961, 0.3921568, 0.1843),3:(0.70588,0.36078, 0.82745),4:(0.52941, 0.831373, 0.262745),5:(0.337254, 0.21176, 0.37647),6:(0.34902, 0.45098, 0.20392),7:(0.32941, 0.21176, 0.145098),8:(0.796078, 0.647059, 0.58431),9:(0.53333, 0.545098, 0.8),10:(0.4941176, 0.82353, 0.545098),11:(0.784314, 0.32941,0.58431),12:(0.3098, 0.415686, 0.41961), 13:(0.6, 0.6, 0.6)}
        self.residues_within_helix={}
        self.universe = topol_object
    def define_amino_acids(self):
        for residue in self.universe.dict_of_plotted_res.keys():
            for aa_type, aa in self.amino_acids.items():
                for amino_acid in aa:
                    if residue[0:3] == amino_acid:
                        self.amino_acid_type[residue]=aa_type
        print "Defining amino acid type..."
    def define_helices(self, helices_text_file, offset):
        with open (helices_text_file, "r") as h:
            lines = h.readlines()
            for line in lines:
                helices[int(line.rsplit(",",3)[0])]=[int(line.rsplit(",",3)[1]),int(line.rsplit(",",3)[2])]

        i=0
        for res in helices:
            i+=1
            while i< len(helices):
                assert helices[res+1][0]-helices[res][1] > 0, "Helices are overlapping in your input file "+helices_text_file+". Please check!"
                break
    
        for helix in helices:
            for res in self.universe.dict_of_plotted_res:
                if int(res[3::])-int(offset)>=helices[helix][0] and int(res[3::])-int(offset)<=helices[helix][1]:
                    self.residues_within_helix[res]=helix
        print "Defining helices..."
        
    def plot_amino_diagramms(self):
        for res in self.amino_acid_type:
            color = [self.colors_amino_acids[self.amino_acid_type[res]],'white']
            plt.figure(figsize=(1.5,1.5))
            #plot a random value (1) at the moment
            ring1,_=plt.pie([1],  radius=1, startangle=90, colors=color, counterclock=False)
            plt.axis('equal')
            plt.setp(ring1, width=1, edgecolor=self.colors_amino_acids[self.amino_acid_type[res]])
            plt.text(0,-0.55,res[0:3]+"\n"+res[3::],ha='center',size=24, fontweight="bold")
            pylab.savefig(str(res[3::])+".svg", dpi=100, transparent=True)
        
        print "Plotting..."
        
    def plot_helices_diagramms(self):
        width=0.20        
        for res in self.universe.closest_atoms:
            if res in self.residues_within_helix:
                color = [self.colors_helices[self.residues_within_helix[res]],'white']
                plt.figure(figsize=(1.5,1.5))
                ring1,_=plt.pie([1],  radius=1-width, startangle=90, colors=color, counterclock=False)
                plt.axis('equal')
                plt.setp(ring1, width=width, edgecolor='white')
                plt.text(0,-0.4,res[0:3]+"\n"+res[3::],ha='center',size=16, fontweight='bold')  
                pylab.savefig(str(res[3::])+".svg", dpi=100, transparent = True)
            else:
                color = [self.colors_helices[13], 'white'] # Matplotlib colors takes more than 1 color, so a second color must be added
                plt.figure(figsize=(1.5,1.5), dpi=100)
                ring1,_=plt.pie([1],  radius=1-width, startangle=90, colors=color, counterclock=False)
                plt.axis('equal')
                plt.setp(ring1, width=width, edgecolor='white')
                plt.text(0,-0.4,res[0:3]+"\n"+res[3::],ha='center',size=16, fontweight='bold')  
                pylab.savefig(str(res[3::])+".svg", dpi=100, transparent=True)
        print "Plotting..."

    def plot_clock_diagramms(self):
        """Uses matplotlib to plot clock diagramms used for data analysis of trajectories, for example, occurance time  over timecourse of simulation"""
        #Plot the residues in clock diagramm fashion
        colors_1=['#1f78b4','white']
        colors_2=['#33a02c','white']
        colors_3=['#6a3d9a','white']

        for res in self.universe.dict_of_plotted_res.keys():
            if len(self.universe.dict_of_plotted_res[res])==4:
                plt.figure(figsize=(1.6, 1.6))
            else:
                plt.figure(figsize=(1.5, 1.5))
            if [sum(1 for x in v if x) for k,v in self.universe.dict_of_plotted_res.items()][0]==2:
                width=0.25
                ring1,_=plt.pie([self.universe.dict_of_plotted_res[res][1],100-self.universe.dict_of_plotted_res[res][1]],  radius=1-width, startangle=90, colors=colors_1, counterclock=False)
            elif [sum(1 for x in v if x) for k,v in self.universe.dict_of_plotted_res.items()][0]==3:
                width=0.25
                ring1,_=plt.pie([self.universe.dict_of_plotted_res[res][1],100-self.universe.dict_of_plotted_res[res][1]],  radius=1-width, startangle=90, colors=colors_1, counterclock=False)
                ring2,_=plt.pie([self.universe.dict_of_plotted_res[res][2],100-self.universe.dict_of_plotted_res[res][2]],  radius=1, startangle=90, colors=colors_2, counterclock=False)
            elif [sum(1 for x in v if x) for k,v in self.universe.dict_of_plotted_res.items()][0]==4:
                width=0.25
                ring1,_=plt.pie([self.universe.dict_of_plotted_res[res][1],100-self.universe.dict_of_plotted_res[res][1]],  radius=1-width, startangle=90, colors=colors_1, counterclock=False)
                ring2,_=plt.pie([self.universe.dict_of_plotted_res[res][2],100-self.universe.dict_of_plotted_res[res][2]],  radius=1, startangle=90, colors=colors_2, counterclock=False)
                ring3,_=plt.pie([self.universe.dict_of_plotted_res[res][3],100-self.universe.dict_of_plotted_res[res][3]],  radius=1+width, startangle=90, colors=colors_3, counterclock=False)
            plt.axis('equal')
            if len(self.universe.dict_of_plotted_res[res])==2:
                plt.setp(ring1, width=width)
                plt.text(0,-0.35,res[0:3]+"\n"+res[3::],ha='center',size=19, fontweight='bold')
            elif len(self.universe.dict_of_plotted_res[res])==3:
                plt.setp(ring1+ring2, width=width)
                plt.text(0,-0.325,res[0:3]+"\n"+res[3::],ha='center',size=15, fontweight='bold')
            elif len(self.universe.dict_of_plotted_res[res])==4:
                plt.setp(ring1+ring2+ring3, width=width)
                plt.text(0,-0.3,res[0:3]+"\n"+res[3::],ha='center',size=13, fontweight='bold')
            pylab.savefig(str(res[3::])+".svg", dpi=100, transparent=True)
            #plt.show()
        print "Plotting..."
