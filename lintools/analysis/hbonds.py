from topol import Topol_Data
import MDAnalysis
import rdkit
from rdkit import Chem
from MDAnalysis.analysis import hbonds
import numpy as np

class HBonds(object):
    def __init__(self,topol_object, topology, trajectory, ligand_name,offset,frame_cutoff,mol2_input = None,pdb_input= None):
        self.HDonorSmarts = Chem.MolFromSmarts('[$([N;!H0;v3]),$([N;!H0;+1;v4]),$([O,S;H1;+0]),$([n;H1;+0])]')
        haccep = "[$([O,S;H1;v2]-[!$(*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),$([N;v3;!$(N-*=!@[O,N,P,S])]),$([nH0,o,s;+0])]"
        self.HAcceptorSmarts = Chem.MolFromSmarts(haccep) 
        self.donors = []
        self.acceptors = []
        self.universe=topol_object
        self.h_bonds = None
        self.mol2_input = mol2_input
        self.pdb_input =pdb_input
        self.ligand = None
        self.hbonds_for_drawing = []
        self.find_donors_and_acceptors_in_ligand()
        self.analyse_hbonds(topology, trajectory, ligand_name,offset,frame_cutoff,3.5)
    def find_donors_and_acceptors_in_ligand(self):
        atom_names=[x.name for x in self.universe.ligand]
        self.ligand = Chem.MolFromMol2File(self.mol2_input,removeHs=False)
        if self.ligand == None:
            print "Exiting. No mol2 file was supplied."
            sys.exit()
        try:
            for atom in self.ligand.GetSubstructMatches(self.HDonorSmarts, uniquify=1):
                self.donors.append(atom_names[atom[0]])
        except AttributeError:
            self.ligand = Chem.MolFromMol2File(self.universe.mol2_file,removeHs=False,sanitize=False)
            self.ligand.UpdatePropertyCache(strict=False)
            for atom in self.ligand.GetSubstructMatches(self.HDonorSmarts, uniquify=1):
                self.donors.append(atom_names[atom[0]])

        for atom in self.ligand.GetSubstructMatches(self.HAcceptorSmarts, uniquify=1):
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
            try:
                md_sim = Topol_Data(topology,None,ligand_name,offset,self.mol2_input)
                h = MDAnalysis.analysis.hbonds.HydrogenBondAnalysis(md_sim.universe,prot_sel[:-3],'(segid '+str(self.universe.ligand.segids[0])+' and resid '+str(self.universe.ligand.resids[0])+')',distance=distance,acceptors=self.acceptors,donors=self.donors)
                h.run()
                h.generate_table()  
                self.h_bonds=h.table
                self.count_hbond_freq(0)

            except ValueError:
                md_sim = Topol_Data(topology,None,ligand_name,offset,self.mol2_input)
                #The curious case of offending residue names that include numbers 
                test = md_sim.universe.select_atoms('(segid '+str(self.universe.ligand.segids[0])+' and resid '+str(self.universe.ligand.resids[0])+')')
                test.resnames = "LIG"
                h = MDAnalysis.analysis.hbonds.HydrogenBondAnalysis(md_sim.universe,prot_sel[:-3],'(segid '+str(self.universe.ligand.segids[0])+' and resid '+str(self.universe.ligand.resids[0])+')',distance=distance,acceptors=self.acceptors,donors=self.donors)
                h.run()
                h.generate_table()  
                self.h_bonds=h.table 
                self.count_hbond_freq(0)

        else:
            try:
                i=0
                for traj in trajectory:
                    md_sim = Topol_Data(topology,trajectory[i],ligand_name, offset,self.mol2_input)
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
                        if self.hbond_frequency[i][bond]>self.universe.frame_count[i]*len(traj)*frame_cutoff/100:
                            self.hbonds_for_drawing.append(bond)
            i+=1

    def count_hbond_freq(self,traj):
        ligand_universe = MDAnalysis.Universe(self.universe.pdb)
        ligand = self.universe.universe.select_atoms('(segid '+str(self.universe.ligand.segids[0])+' and resid '+str(self.universe.ligand.resids[0])+')')
        for i in range(np.prod(self.h_bonds.shape)):
            if self.h_bonds[i][3]==ligand.resnames[0]:
                atomname = self.h_bonds[i][5]
            else:
                atomname = self.h_bonds[i][8]
            if atomname.startswith("O",0):
                lig_atom=atomname
            else:
                for atom in ligand_universe.atoms:
                    if atomname == atom.name:
                        atom_id = atom.id
                rdkit_atom = self.ligand.GetAtomWithIdx(atom_id)
                for neigh in rdkit_atom.GetNeighbors():
                    neigh_atom_id = neigh.GetIdx()
                lig_atom = ligand_universe.atoms[neigh_atom_id].name
                #find the idx of this atomname 
            if self.h_bonds[i][3]==ligand.resnames[0]:
                results_tuple = lig_atom,self.h_bonds[i][6]+str(self.h_bonds[i][7])
            else:
                results_tuple = lig_atom,self.h_bonds[i][3]+str(self.h_bonds[i][4])
            if results_tuple not in self.hbond_frequency.keys():
                self.hbond_frequency[traj][results_tuple]=1
            else:
                self.hbond_frequency[traj][results_tuple]=int(self.hbond_frequency[traj][results_tuple])+1


