import mdtraj
from timeit import default_timer as timer
from collections import namedtuple, defaultdict


class SASA(object):
	"""This module analyses the solvent accessible surface area (SASA) of the ligand
	using a module from MDTraj. It calculates SASA for every atom during the simulation
	and returns all the values as well as averaged values.
	Takes:
		* topology data object * - information about system (lintools.Data object)
		* trajectory * - list of trajectories
	Output:
		* atom_sasa * - ligand atom SASA
		* total_sasa * - averaged SASA of the ligand atoms over the time of the simulation
		"""
	def __init__(self,topology_data_object,trajectory):
		self.topology_data = topology_data_object
		self.trajectory = trajectory

		self.atom_sasa = {}
		self.analyse_ligand_sasa()
	def analyse_ligand_sasa(self):
		"""Analysis of ligand SASA."""
		i=0
		start = timer()
		if self.trajectory == []:
			self.trajectory = [self.topology_data.universe.filename]
		try:
			for traj in self.trajectory:
				new_traj = mdtraj.load(traj,top=self.topology_data.universe.filename)
				#Analyse only non-H ligand
				ligand_slice = new_traj.atom_slice(atom_indices=self.topology_data.universe.ligand_noH.ids)

				self.sasa = mdtraj.shrake_rupley(ligand_slice)
				self.atom_sasa[i]=self.assign_per_atom_sasa()
				i+=1
			self.total_sasa = self.get_total_per_atom_sasa()
		except KeyError as e:
			print "WARNING: SASA analysis cannot be performed due to incorrect atom names in"
			print "the topology ", e

		print "SASA: "+str(timer()-start)


	def assign_per_atom_sasa(self):
		"""Make a dictionary with SASA assigned to each ligand atom, stored as list of SASA values over
		the simulation time."""
		atom_names= [atom.name for atom in self.topology_data.universe.ligand_noH.atoms]
		sasa_dict = {}
		for atom in range(0,self.topology_data.universe.ligand_noH.n_atoms):
			sasa_dict[atom_names[atom]]=[self.sasa[i][atom] for i in range(len(self.sasa))]
		return sasa_dict

	def get_total_per_atom_sasa(self):
		"""Return average SASA of the atoms."""
		total_sasa = defaultdict(int)
		for traj in range(len(self.atom_sasa)):
			for atom in self.atom_sasa[traj]:
				total_sasa[atom]+=float(sum((self.atom_sasa[traj][atom])))/len(self.atom_sasa[traj][atom])
		for atom in total_sasa:
			total_sasa[atom]=float(total_sasa[atom])/len(self.atom_sasa)
		return total_sasa
