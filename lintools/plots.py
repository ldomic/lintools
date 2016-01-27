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
        self.amino_acids = {"acidic":["ASP","GLU"], "basic":["LYS","ARG"], "aromatic":["PHE","TYR","TRP"],"polar":["SER","THR","ASN","GLN","CYS","HIS"],"hydrophobic":["ALA","VAL","ILE","LEU","MET","GLY","PRO"]}
        self.colors_amino_acids = {"acidic":"#D9774B", "basic":"#889DCC", "aromatic":"#9FC74A", "polar":"#D06AC1","hydrophobic":"#6AC297"}
        self.amino_acid_type={}
        self.colors_domains ={1:"#78D035",2:"#BE4F24",3:"#7676C2",4:"#7CC2C0",5:"#42291B",6:"#BB3B83",7:"#486128",8:"#6FCD7A",9:"#A43FC8",10:"#C09584",11:"#C5AD3D",12:"#43264D", 13:"#A9A9A9"}
        self.residues_within_domain={}
        self.plotted_domains = []
        self.universe = topol_object
    def define_amino_acids(self):
        for residue in self.universe.dict_of_plotted_res.keys():
            for aa_type, aa in self.amino_acids.items():
                for amino_acid in aa:
                    if residue[0:3] == amino_acid:
                        self.amino_acid_type[residue]=aa_type
        print "Defining amino acid type..."
    def define_domains(self, domain_file, offset):
        domains = {}
        with open (domain_file, "r") as h:
            lines = h.readlines()
            for line in lines:
                if len(line.rsplit(";",5))==3:
                    domains[int(line.rsplit(";",5)[0])]=([],line.rsplit(";",5)[2])
                if len(line.rsplit(";",5))==4:
                    domains[int(line.rsplit(";",5)[0])]=([],line.rsplit(";",5)[2],line.rsplit(";",5)[3][:7])
                if len(line.rsplit(";",5))==5:
                    domains[int(line.rsplit(";",5)[0])]=([],line.rsplit(";",5)[2],line.rsplit(";",5)[3][:7],line.rsplit(";",5)[4])
                for i in range(len(line.rsplit(";",5)[1].rsplit(",",50))):
                    if len(line.rsplit(";",5)[1].rsplit(",",50)[i].rsplit("-",2))>1:
                        for num in range(int(line.rsplit(";",5)[1].rsplit(",",50)[i].rsplit("-",2)[0]),int(line.rsplit(";",5)[1].rsplit(",",50)[i].rsplit("-",2)[1])+1):
                            domains[int(line.rsplit(";",5)[0])][0].append(num)
                    else:
                        domains[int(line.rsplit(";",5)[0])][0].append(int(line.rsplit(";",5)[1].rsplit(",",50)[i]))
        for domain in domains:
            for res in self.universe.dict_of_plotted_res:
                if int(res[3::]) in domains[domain][0] and len(domains[domain])==2:
                    self.residues_within_domain[res]=[domain,domains[domain][1],self.colors_domains[domain],"N"]
                if int(res[3::]) in domains[domain][0] and len(domains[domain])==3:
                    self.residues_within_domain[res]=[domain,domains[domain][1],domains[domain][2],"N"]
                if int(res[3::]) in domains[domain][0] and len(domains[domain])==4:
                    #extra checks for dashes
                    if str(domains[domain][2])=="None":
                        print res
                        if res not in self.residues_within_domain.keys():
                            self.residues_within_domain[res]=[domain, domains[domain][1],"#A9A9A9",domains[domain][3]]
                        self.residues_within_domain[res][3]=[domains[domain][3]]
                        print self.residues_within_domain[res]
                        domain_description=[domain, domains[domain][1],"#A9A9A9",domains[domain][3]]
                        if domain_description not in self.plotted_domains:
                            #domain_description=[domain, domains[domain][1],"#A9A9A9",domains[domain][3]]
                            self.plotted_domains.append(domain_description)
                    else:
                        self.residues_within_domain[res]=[domain,domains[domain][1],domains[domain][2], domains[domain][3]]


        for res in self.universe.dict_of_plotted_res:
            if res not in self.residues_within_domain.keys():
                self.residues_within_domain[res]=[0,"No specified domain", "#A9A9A9","N"]
        for res in self.residues_within_domain:
            if self.residues_within_domain[res] not in self.plotted_domains:
                self.plotted_domains.append(self.residues_within_domain[res])

        print "Defining domains..."
        
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
        
    def plot_domains_diagramms(self):
        width=0.20        
        for res in self.universe.closest_atoms:
            if res in self.residues_within_domain:
                color = [self.residues_within_domain[res][2],'white']
                plt.figure(figsize=(1.5,1.5))
                ring1,_=plt.pie([1],  radius=1-width, startangle=90, colors=color, counterclock=False)
                plt.axis('equal')
                plt.setp(ring1, width=width, edgecolor='white')
                plt.text(0,-0.35,res[0:3]+"\n"+res[3::],ha='center',size=20, fontweight='bold')  
                pylab.savefig(str(res[3::])+".svg", dpi=100, transparent = True)
            else:
                color = [self.colors_domains[13], 'white'] # Matplotlib colors takes more than 1 color, so a second color must be added
                plt.figure(figsize=(1.5,1.5), dpi=100)
                ring1,_=plt.pie([1],  radius=1-width, startangle=90, colors=color, counterclock=False)
                plt.axis('equal')
                plt.setp(ring1, width=width, edgecolor='white')
                plt.text(0,-0.35,res[0:3]+"\n"+res[3::],ha='center',size=20, fontweight='bold')  
                pylab.savefig(str(res[3::])+".svg", dpi=100, transparent=True)
        print "Plotting..."

    def plot_clock_diagramms_old(self):
        """Uses matplotlib to plot clock diagramms used for data analysis of trajectories, for example, occurrence time  over timecourse of simulation"""
        #Plot the residues in clock diagramm fashion
        colors_1=['#1f78b4','white']
        #colors_1=['#5C0016','white']
        colors_2=['#33a02c','white']
        colors_3=['#6a3d9a','white']

        for res in self.universe.dict_of_plotted_res.keys():
            if len(self.universe.dict_of_plotted_res[res])==4:
                plt.figure(figsize=(1.6, 1.6))
            else:
                plt.figure(figsize=(1.5, 1.5))
            if [sum(1 for x in v if x) for k,v in self.universe.dict_of_plotted_res.items()][0]==2:
                width=0.25
                ring1,_=plt.pie([self.universe.dict_of_plotted_res[res][1],self.universe.frame_count-self.universe.dict_of_plotted_res[res][1]],  radius=1-width, startangle=90, colors=colors_1, counterclock=False)
            elif [sum(1 for x in v if x) for k,v in self.universe.dict_of_plotted_res.items()][0]==3:
                width=0.25
                ring1,_=plt.pie([self.universe.dict_of_plotted_res[res][1],self.universe.frame_count-self.universe.dict_of_plotted_res[res][1]],  radius=1-width, startangle=90, colors=colors_1, counterclock=False)
                ring2,_=plt.pie([self.universe.dict_of_plotted_res[res][2],self.universe.frame_count-self.universe.dict_of_plotted_res[res][2]],  radius=1, startangle=90, colors=colors_2, counterclock=False)
            elif [sum(1 for x in v if x) for k,v in self.universe.dict_of_plotted_res.items()][0]==4:
                width=0.25
                ring1,_=plt.pie([self.universe.dict_of_plotted_res[res][1],self.universe.frame_count-self.universe.dict_of_plotted_res[res][1]],  radius=1-width, startangle=90, colors=colors_1, counterclock=False)
                ring2,_=plt.pie([self.universe.dict_of_plotted_res[res][2],self.universe.frame_count-self.universe.dict_of_plotted_res[res][2]],  radius=1, startangle=90, colors=colors_2, counterclock=False)
                ring3,_=plt.pie([self.universe.dict_of_plotted_res[res][3],self.universe.frame_count-self.universe.dict_of_plotted_res[res][3]],  radius=1+width, startangle=90, colors=colors_3, counterclock=False)
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
    def plot_clock_diagramms(self):
        """Uses matplotlib to plot clock diagramms used for data analysis of trajectories, for example, occurrence time  over timecourse of simulation"""
        colors = [['#00441b','white'], ['#1b7837','white'],['#5aae61','white'], ['#9970ab','white'], ['#762a83','white'],['#40004b',"white"]]
        for res in self.universe.dict_of_plotted_res.keys():
            plt.figure(figsize=(1.55, 1.55))
            ring_number=[sum(1 for x in v if x) for k,v in self.universe.dict_of_plotted_res.items()][0]
            width=0.75/ring_number
            rings=[]
            for ring in range(1,ring_number):
                ring,_=plt.pie([self.universe.dict_of_plotted_res[res][ring],self.universe.frame_count-self.universe.dict_of_plotted_res[res][ring]],  radius=0.65+width*ring, startangle=90, colors=colors[ring-1], counterclock=False)
                rings=rings+ring
            plt.setp(rings, width=width)
            plt.text(0,-0.3,res[0:3]+"\n"+res[3::],ha='center',size=13, fontweight='bold')
            pylab.savefig(str(res[3::])+".svg", dpi=100, transparent=True)

