from topol import Topol_Data
import MDAnalysis
import rdkit
from rdkit import Chem
from MDAnalysis.analysis import hbonds
import numpy as np

class HBonds(object):
    def __init__(self,topol_object, mol2_file, topology, trajectory, ligand_name,offset,frame_cutoff):
        self.HDonorSmarts = Chem.MolFromSmarts('[$([N;!H0;v3]),$([N;!H0;+1;v4]),$([O,S;H1;+0]),$([n;H1;+0])]')
        haccep = "[$([O,S;H1;v2]-[!$(*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),$([N;v3;!$(N-*=!@[O,N,P,S])]),$([nH0,o,s;+0])]"
        self.HAcceptorSmarts = Chem.MolFromSmarts(haccep) 
        self.donors = []
        self.acceptors = []
        self.universe=topol_object
        self.h_bonds = None
        self.hbonds_for_drawing = []
        self.mol2_file = mol2_file
        self.find_donors_and_acceptors_in_ligand()
        self.analyse_hbonds(topology, trajectory, ligand_name,offset,frame_cutoff)
    def find_donors_and_acceptors_in_ligand(self):
        atom_names=[x.name for x in self.universe.ligand]
        ligand = Chem.MolFromMol2File(self.mol2_file,removeHs=False)
        for atom in ligand.GetSubstructMatches(self.HDonorSmarts, uniquify=1):
            self.donors.append(atom_names[atom[0]])
        for atom in ligand.GetSubstructMatches(self.HAcceptorSmarts, uniquify=1):
             self.acceptors.append(atom_names[atom[0]])
    def analyse_hbonds(self,topology, trajectory, ligand_name,offset,frame_cutoff):
        prot_sel = "protein and "
        for res in self.universe.dict_of_plotted_res.values():
            prot_sel=prot_sel+"resid "+str(res[0])+" or "
        if trajectory is None:
            md_sim = Topol_Data(topology,None,ligand_name,offset)
            h = MDAnalysis.analysis.hbonds.HydrogenBondAnalysis(md_sim.universe,prot_sel[:-3],'(segid '+str(md_sim.ligand.segids[0])+' and resid '+str(md_sim.ligand.resids[0])+')',acceptors=self.acceptors,donors=self.donors)
            h.run()
            h.generate_table()  
            self.h_bonds=h.table  
        else:
            md_sim = Topol_Data(topology,None,ligand_name,offset)
            h = MDAnalysis.analysis.hbonds.HydrogenBondAnalysis(md_sim.universe,prot_sel[:-3],'(segid '+str(md_sim.ligand.segids[0])+' and resid '+str(md_sim.ligand.resids[0])+')',acceptors=self.acceptors,donors=self.donors)
            h.run()
            h.generate_table()  
            self.h_bonds=h.table  
            i=0
            for traj in trajectory:
                i+=1
                md_sim = Topol_Data(topology,traj,ligand_name, offset)
                h = MDAnalysis.analysis.hbonds.HydrogenBondAnalysis(md_sim.universe,prot_sel[:-3],'(segid '+str(ligand_name.segids[0])+' and resid '+str(ligand_name.resids[0])+')',acceptors=self.acceptors,donors=self.donors)
                h.run()
                h.generate_table()  
                print h.table.shape
                self.h_bonds=np.hstack((self.h_bonds,h.table)) 
        self.distance = h.distance
        #h = MDAnalysis.analysis.hbonds.HydrogenBondAnalysis(self.universe.universe,"protein",'(segid '+str(self.universe.ligand.segids[0])+' and resid '+str(self.universe.ligand.resids[0])+')',acceptors=self.acceptors,donors=self.donors)
        ligand_from_mol2 = MDAnalysis.Universe(self.mol2_file) 
        self.hbond_frequency = {}
        for i in range(np.prod(self.h_bonds.shape)):
            if self.h_bonds[i][3]==self.universe.ligand.resnames[0]:
                atomname = self.h_bonds[i][5]
            else:
                atomname = self.h_bonds[i][8]
            if atomname.startswith("O",0):
                lig_atom=atomname
            else:
                for x in range(len(ligand_from_mol2.atoms.bonds.bondlist)):
                    if atomname==ligand_from_mol2.atoms.bonds.bondlist[x][0].name:
                        lig_atom =  ligand_from_mol2.atoms.bonds.bondlist[x][1].name
                    if atomname==ligand_from_mol2.atoms.bonds.bondlist[x][1].name:
                        lig_atom = ligand_from_mol2.atoms.bonds.bondlist[x][0].name
            if self.h_bonds[i][3]==self.universe.ligand.resnames[0]:
                results_tuple = lig_atom,self.h_bonds[i][6]+str(self.h_bonds[i][7])
            else:
                results_tuple = lig_atom,self.h_bonds[i][3]+str(self.h_bonds[i][4])
            if results_tuple not in self.hbond_frequency.keys():
                self.hbond_frequency[results_tuple]=1
            else:
                self.hbond_frequency[results_tuple]=int(self.hbond_frequency[results_tuple])+1
        for bond in self.hbond_frequency:
            if trajectory!=None:
                if self.hbond_frequency[bond]>self.universe.frame_count*len(trajectory)*frame_cutoff/100:
                    self.hbonds_for_drawing.append(bond)
            else:
                if self.hbond_frequency[bond]>self.universe.frame_count*frame_cutoff/100:
                    self.hbonds_for_drawing.append(bond)