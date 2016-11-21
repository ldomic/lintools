from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import  rdPartialCharges
from collections import namedtuple

class LigDescr(object):
	def __init__(self,topology_data_object):
		self.topology_data = topology_data_object
		mol = Chem.RemoveHs(self.topology_data.mol2)
		self.calculate_descriptors(mol)
	def calculate_descriptors(self,mol):
		#make dictionary
		self.ligand_atoms = {index:{"name":x.name} for index,x in enumerate(self.topology_data.universe.ligand_noH.atoms)}

		#Calculate logP and MR
		contribs = self.calculate_logP(mol)

		#Calculate Gasteiger charges
		self.calculate_Gasteiger_charges(mol)

		#Calculate formal charges
		fcharges = self.calculate_formal_charge(mol)

		for atom in self.ligand_atoms.keys():
			self.ligand_atoms[atom]["logP"]=contribs[atom][0]
			self.ligand_atoms[atom]["MR"]=contribs[atom][1]
			self.ligand_atoms[atom]["Gasteiger_ch"]=mol.GetAtomWithIdx(atom).GetProp("_GasteigerCharge")
			self.ligand_atoms[atom]["Formal charges"]=fcharges[atom]

		#Determine rotatable bonds
		self.rot_bonds=self.get_rotatable_bonds(mol)


	def calculate_logP(self,mol):
		"""Calculates Crippen contributions, i.e. logP of ligand molecule.
		Takes:
			* mol * - mol2 file in rdkit environment
		Returns:
			* contribs * - tuple of Wildman-Crippen logP, MR (molar refractivity - 
				measure of the volume occupied by a molecule of the substance) values
		"""
		contribs = rdMolDescriptors._CalcCrippenContribs(mol)
		return contribs
	def get_rotatable_bonds(self,mol):
		"""Determines rotatable bonds in a ligand molecule
		Takes:
			* mol * - mol2 file in rdkit environment
		Output:
			* bonds * - tuples of atom ids 
		"""
		RotatableBondSmarts = Chem.MolFromSmarts('[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]') 
		bonds = mol.GetSubstructMatches(RotatableBondSmarts,uniquify=1)
		return bonds
	def calculate_Gasteiger_charges(self,mol):
		rdPartialCharges.ComputeGasteigerCharges(mol)
	def calculate_formal_charge(self,mol):
		formal_charges = []
		for atom in mol.GetAtoms():
			formal_charges.append(atom.GetFormalCharge())
		return formal_charges

