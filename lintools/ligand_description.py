from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import  rdPartialCharges
from collections import namedtuple

class LigDescr(object):
	__version__ = "09.2016"
	"""This module analyses various properties of the ligand molecule such as charges, logP and MR. 
	It also detects rotatable bonds.
	Takes:
		* topology data object * - information about system (lintools.Data object)
	Output:
		* ligand_atoms * - dictionary with information about ligand properties
		* rot_bonds * - list of tuples with atom indices
		"""
	def __init__(self,topology_data_object,rmsf_object,sasa_object):
		self.topology_data = topology_data_object
		self.sasa = sasa_object
		self.rmsf = rmsf_object
		mol = Chem.RemoveHs(self.topology_data.mol2)
		self.calculate_descriptors(mol)
	def calculate_descriptors(self,mol):
		"""Calculates descriptors such as logP, charges and MR and saves that in a dictionary."""
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
			if self.sasa!=None:
				self.ligand_atoms[atom]["SASA"]=self.sasa.total_sasa[self.ligand_atoms[atom]["name"]]
			if self.rmsf!=None:
				self.ligand_atoms[atom]["RMSF"]=self.rmsf.ligand_rmsf[self.ligand_atoms[atom]["name"]]

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
		"""Calculates Gasteiger partial charges saved as prop under the name (_GasteigerCharge)
		"""
		rdPartialCharges.ComputeGasteigerCharges(mol)
	def calculate_formal_charge(self,mol):
		"""Calculates formal charge for each atom.
		Takes:
			* mol * - mol2 file in rdkit environment
		Output:
			* formal_charges * - list of charges
		"""
		formal_charges = []
		for atom in mol.GetAtoms():
			formal_charges.append(atom.GetFormalCharge())
		return formal_charges

