import MDAnalysis
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.pyplot import cm 

class Dimensionality_Reduction(object):
	def __init__(self, topology, ligand_name,cutoff=3.5):
		self.topology = topology
		self.cutoff = cutoff
		self.ligand = ligand_name
		self.make_coords()
	def reduce_dimensions(self, a):
		mean_x = np.mean(a[0,:])
		mean_y = np.mean(a[1,:])
		mean_z = np.mean(a[2,:])

		mean_vector = np.array([[mean_x],[mean_y],[mean_z]])
		#calculate d dimensional mean vector
		scatter_matrix = np.zeros((3,3))
		for i in range(a.shape[1]):
		    scatter_matrix += (a[:,i].reshape(3,1) - mean_vector).dot((a[:,i].reshape(3,1) - mean_vector).T)
		    
		eig_val_sc, eig_vec_sc = np.linalg.eig(scatter_matrix)

		# Make a list of (eigenvalue, eigenvector) tuples
		eig_pairs = [(np.abs(eig_val_sc[i]), eig_vec_sc[:,i]) for i in range(len(eig_val_sc))]

		# Sort the (eigenvalue, eigenvector) tuples from high to low
		eig_pairs.sort()
		eig_pairs.reverse()

		matrix_w = np.hstack((eig_pairs[0][1].reshape(3,1), eig_pairs[1][1].reshape(3,1)))
		transformed = matrix_w.T.dot(a)
		return transformed

	def make_coords(self):
		my_univ = MDAnalysis.Universe(self.topology)
	    	lig = my_univ.select_atoms("resname "+self.ligand+" and not name H*")
    		closest_res = my_univ.select_atoms("protein and around "+str(self.cutoff)+" resname "+self.ligand)
    		closest_res_ca = closest_res.residues.atoms["CA"]
    		res_pos = np.transpose(closest_res_ca.positions)
    		lig_pos = np.transpose(lig.positions)
    		all_pos = np.hstack((res_pos, lig_pos))
    		all_new_pos = self.reduce_dimensions(all_pos)
    		i=0
    		color=iter(cm.rainbow(np.linspace(0,1,len(closest_res_ca))))
		for res in closest_res_ca:
		    c=next(color)
		    plt.plot(all_new_pos[0,i:i+1], all_new_pos[1,i:i+1], 'o', markersize=12, c=c, alpha=0.5)
		    plt.annotate(res.resname+str(res.resid),xy=(all_new_pos[0,i:i+1], all_new_pos[1,i:i+1]))
		    i+=1
		plt.plot(all_new_pos[0,len(closest_res_ca):], all_new_pos[1,len(closest_res_ca):], 'o', markersize=7, color='blue', alpha=0.5, label=self.ligand)
		plt.xlabel('x_values')
		plt.ylabel('y_values')
		plt.legend(loc=9,bbox_to_anchor=(0.5, -0.1))

		plt.show()

