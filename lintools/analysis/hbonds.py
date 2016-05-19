from topol import Topol_Data
import MDAnalysis
import rdkit
from rdkit import Chem
from MDAnalysis.analysis import hbonds
import numpy as np

class HBonds(object):
    def __init__(self,topol_object, topology, trajectory, ligand_name,offset,frame_cutoff,tests=False):
        self.HDonorSmarts = Chem.MolFromSmarts('[$([N;!H0;v3]),$([N;!H0;+1;v4]),$([O,S;H1;+0]),$([n;H1;+0])]')
        haccep = "[$([O,S;H1;v2]-[!$(*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),$([N;v3;!$(N-*=!@[O,N,P,S])]),$([nH0,o,s;+0])]"
        self.HAcceptorSmarts = Chem.MolFromSmarts(haccep) 
        self.donors = []
        self.acceptors = []
        self.universe=topol_object
        self.h_bonds = None
        self.hbonds_for_drawing = []
        if tests==False:
            self.find_donors_and_acceptors_in_ligand()
            self.analyse_hbonds(topology, trajectory, ligand_name,offset,frame_cutoff,3.5)
    def find_donors_and_acceptors_in_ligand(self):
        atom_names=[x.name for x in self.universe.ligand]
        for atom in self.universe.mol2.GetSubstructMatches(self.HDonorSmarts, uniquify=1):
            self.donors.append(atom_names[atom[0]])
        for atom in self.universe.mol2.GetSubstructMatches(self.HAcceptorSmarts, uniquify=1):
            self.acceptors.append(atom_names[atom[0]])
    def analyse_hbonds(self,topology, trajectory, ligand_name,offset,frame_cutoff,distance=3.5):
        self.hbond_frequency={}
        if trajectory!=None:
            i=0
            for traj in trajectory:
                self.hbond_frequency[i]={}
                i+=1
        else:
            self.hbond_frequency[0]={}
        prot_sel = "protein and "
        for res in self.universe.dict_of_plotted_res.values():
            prot_sel=prot_sel+"resid "+str(res[0])+" or "
        if trajectory is None:
            md_sim = Topol_Data(topology,None,ligand_name,offset)
            ligand = md_sim.universe.select_atoms('(segid '+str(self.universe.ligand.segids[0])+' and resid '+str(self.universe.ligand.resids[0])+')')
            ligand.resnames = "LIG"
            ligand.resname = "LIG"
            h = MDAnalysis.analysis.hbonds.HydrogenBondAnalysis(md_sim.universe,prot_sel[:-3],'(segid '+str(self.universe.ligand.segids[0])+' and resid '+str(self.universe.ligand.resids[0])+')',distance=distance,acceptors=self.acceptors,donors=self.donors)
            h.run()
            h.generate_table()  
            self.h_bonds=h.table
            self.count_hbond_freq(0)
        else:
            try:
                i=0
                for traj in trajectory:
                    md_sim = Topol_Data(topology,trajectory[i],ligand_name, offset)
                    h = MDAnalysis.analysis.hbonds.HydrogenBondAnalysis(md_sim.universe,prot_sel[:-3],'(segid '+str(self.universe.ligand.segids[0])+' and resid '+str(self.universe.ligand.resids[0])+')',distance=distance,acceptors=self.acceptors,donors=self.donors)
                    h.run()
                    h.generate_table()  
                    self.h_bonds=h.table
                    self.count_hbond_freq(i)
                    i+=1
            except AttributeError:
                i=0
                for traj in trajectory:
                    md_sim = MDAnalysis.Universe(topology,trajectory[i])
                    ligand = md_sim.universe.select_atoms('(segid '+str(self.universe.ligand.segids[0])+' and resid '+str(self.universe.ligand.resids[0])+')')
                    ligand.resnames = "LIG"
                    ligand.resname = "LIG"
                    h = MDAnalysis.analysis.hbonds.HydrogenBondAnalysis(md_sim.universe,"protein",'(segid '+str(self.universe.ligand.segids[0])+' and resid '+str(self.universe.ligand.resids[0])+')',distance=distance, acceptors=self.acceptors,donors=self.donors)
                    h.run()
                    h.generate_table()  
                    self.h_bonds=h.table
                    self.count_hbond_freq(i)
                    i+=1
            except ValueError:
                i=0
                for traj in trajectory:
                    md_sim = MDAnalysis.Universe(topology,trajectory[i])
                    ligand = md_sim.universe.select_atoms('(segid '+str(self.universe.ligand.segids[0])+' and resid '+str(self.universe.ligand.resids[0])+')')
                    ligand.resnames = "LIG"
                    ligand.resname = "LIG"
                    h = MDAnalysis.analysis.hbonds.HydrogenBondAnalysis(md_sim.universe,"protein",'(segid '+str(self.universe.ligand.segids[0])+' and resid '+str(self.universe.ligand.resids[0])+')', distance=distance,acceptors=self.acceptors,donors=self.donors)
                    h.run()
                    h.generate_table()  
                    self.h_bonds=h.table
                    self.count_hbond_freq(i)
                    i+=1
        self.distance = h.distance
        if trajectory == None:
            for traj in self.hbond_frequency.values():
                for bond in traj:
                    self.hbonds_for_drawing.append(bond)
        else:
            i=0
            for traj in trajectory:
                for the_traj in self.hbond_frequency.values():
                    for bond in the_traj:
                        if self.hbond_frequency[i][bond]>self.universe.frame_count[i]*int(frame_cutoff)/100 and bond not in self.hbonds_for_drawing:
                            self.hbonds_for_drawing.append(bond)
            i+=1

    def count_hbond_freq(self,traj):
        #Has been upgraded to work with changes in MDAnalysis 0.15.0
        ligand_universe = MDAnalysis.Universe(self.universe.pdb)
        ligand = self.universe.universe.select_atoms('(segid '+str(self.universe.ligand.segids[0])+' and resid '+str(self.universe.ligand.resids[0])+')')
        ligand.resnames = "LIG"   
        ligand.resname = "LIG"
        for i in range(np.prod(self.h_bonds.shape)):
            print self.h_bonds
            if self.h_bonds[i][5]==ligand.resnames[0]:
                atomname = self.h_bonds[i][7]
            else:
                atomname = self.h_bonds[i][10]
            if atomname.startswith("O",0) or atomname.startswith("N",0):
                lig_atom=atomname
            else:
                for atom in ligand_universe.atoms:
                    if atomname == atom.name:
                        print atomname, atom.id
                        atom_id = int(atom.id)-1
                rdkit_atom = self.universe.mol2.GetAtomWithIdx(atom_id)
                for neigh in rdkit_atom.GetNeighbors():
                    neigh_atom_id = neigh.GetIdx()
                lig_atom = ligand_universe.atoms[neigh_atom_id].name
                #find the idx of this atomname 
            #check whether the hydrogen bond is formed with side chain or backbone of residue
            if self.h_bonds[i][5]==ligand.resnames[0]:
                if self.h_bonds[i][10]=="O" or self.h_bonds[i][10]=="N" or self.h_bonds[i][7]=="H":
                    results_tuple = lig_atom,self.h_bonds[i][8]+str(self.h_bonds[i][9]),"backbone"
                else:
                    results_tuple = lig_atom,self.h_bonds[i][8]+str(self.h_bonds[i][9]),"sidechain"
            else:
                if self.h_bonds[i][5]=="O" or self.h_bonds[i][7]=="N" or self.h_bonds[i][7]=="H":
                    results_tuple = lig_atom,self.h_bonds[i][5]+str(self.h_bonds[i][6]),"backbone"
                else:
                    results_tuple = lig_atom,self.h_bonds[i][5]+str(self.h_bonds[i][6]),"sidechain"
            if results_tuple not in self.hbond_frequency[traj].keys():
                self.hbond_frequency[traj][results_tuple]=1
            else:
                self.hbond_frequency[traj][results_tuple]=int(self.hbond_frequency[traj][results_tuple])+1
            

