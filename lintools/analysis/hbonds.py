from topol import Topol_Data
import MDAnalysis
import rdkit
from rdkit import Chem
from MDAnalysis.analysis import hbonds
import numpy as np

class HBonds(object):
    def __init__(self,topol_object,ligand_name):
        self.HDonorSmarts = Chem.MolFromSmarts('[$([N;!H0;v3]),$([N;!H0;+1;v4]),$([O,S;H1;+0]),$([n;H1;+0])]')
        haccep = "[$([O,S;H1;v2]-[!$(*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),$([N;v3;!$(N-*=!@[O,N,P,S])]),$([nH0,o,s;+0])]"
        self.HAcceptorSmarts = Chem.MolFromSmarts(haccep) 
        self.donors = []
        self.acceptors = []
        self.h_bonds = None
        self.hbonds_for_drawing = []
        self.universe = topol_object
        self.find_donors_and_acceptors_in_ligand(ligand_name)
        self.analyse_hbonds(ligand_name)
    def find_donors_and_acceptors_in_ligand(self,ligand_name):
        atom_names=[x.name for x in self.universe.ligand]
        
        ligand = Chem.MolFromMol2File(ligand_name+".mol2",removeHs=False)
        for atom in ligand.GetSubstructMatches(self.HDonorSmarts, uniquify=1):
            self.donors.append(atom_names[atom[0]])
        for atom in ligand.GetSubstructMatches(self.HAcceptorSmarts, uniquify=1):
             self.acceptors.append(atom_names[atom[0]])
    def analyse_hbonds(self, ligand_name):
        h = MDAnalysis.analysis.hbonds.HydrogenBondAnalysis(self.universe.universe,"protein","resname "+ligand_name,acceptors=self.acceptors,donors=self.donors)
        h.run()
        h.generate_table()
        self.h_bonds = h.table
        for i in range(np.prod(self.h_bonds.shape)):
            atomname = self.h_bonds[i][5]
            for x in range(len(self.universe.ligand.atoms.bonds.bondlist)):
                if atomname==self.universe.ligand.atoms.bonds.bondlist[x][0].name:
                    lig_atom =  self.universe.ligand.atoms.bonds.bondlist[x][1].name
                if atomname==self.universe.ligand.atoms.bonds.bondlist[x][1].name:
                    lig_atom = self.universe.ligand.atoms.bonds.bondlist[x][0].name
            results_tuple = lig_atom,self.h_bonds[i][6]+str(self.h_bonds[i][7])
            self.hbonds_for_drawing.append(results_tuple)