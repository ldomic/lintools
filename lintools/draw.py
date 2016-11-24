import matplotlib 
import matplotlib.pyplot as plt
from matplotlib.transforms import IdentityTransform,Affine2D

from matplotlib import pylab
from shapely import geometry
from shapely.ops import cascaded_union
import rdkit
from rdkit import Chem
import fileinput
import sys


class Draw(object):
	def __init__(self,topology_data_object,molecule_object,hbond_object,pi_contacts_object,saltbridges_object,lig_descr_object):
		self.topology_data = topology_data_object
		self.molecule = molecule_object
		self.hbonds = hbond_object
		self.pi_contacts = pi_contacts_object
		self.saltbridges = saltbridges_object
		self.lig_descr = lig_descr_object
		self.cloud=""
		self.draw_hydrogen_bonds()
		self.draw_salt_bridges()
		self.draw_pi_contacts()
		self.add_smiles_id()
		self.draw_molecule(None,None,None,None)
	def draw_hydrogen_bonds(self,color="black"):
		"""For each bond that has been determined to be important, a line gets drawn.
		"""
		self.draw_hbonds ="<g transform='translate("+str((self.molecule.x_dim-self.molecule.molsize1)/2)+","+str((self.molecule.y_dim-self.molecule.molsize2)/2)+")'>'"

		if self.hbonds!=None:
			for bond in self.hbonds.hbonds_for_drawing:
				atom = self.topology_data.universe.atoms[bond[0]-1] #zero-based index vs one-based index
				residue = (atom.resname, str(atom.resid), atom.segid)
				if bond[2] in ["N","O","H"]:
					#backbone interactions
					self.draw_hbonds=self.draw_hbonds+"<line x1='"+str(int(self.molecule.nearest_points_coords[residue][0]))+"' y1='"+str(int(self.molecule.nearest_points_coords[residue][1]))+"' x2='"+str(float(self.molecule.ligand_atom_coords_from_diagr[bond[1]][0]))+"' y2='"+str(float(self.molecule.ligand_atom_coords_from_diagr[bond[1]][1]))+"' style='stroke:"+color+";stroke-width:4' />"
				else:    
					#sidechain interactions
					self.draw_hbonds=self.draw_hbonds+"<line stroke-dasharray='5,5'  x1='"+str(int(self.molecule.nearest_points_coords[residue][0]))+"' y1='"+str(int(self.molecule.nearest_points_coords[residue][1]))+"' x2='"+str(float(self.molecule.ligand_atom_coords_from_diagr[bond[1]][0]))+"' y2='"+str(float(self.molecule.ligand_atom_coords_from_diagr[bond[1]][1]))+"' style='stroke:"+color+";stroke-width:4' />"
		self.draw_hbonds+="</g>"

	def draw_pi_contacts(self):
		if self.pi_contacts!=None:
			self.draw_pi_lines="<g transform='translate("+str((self.molecule.x_dim-self.molecule.molsize1)/2)+","+str((self.molecule.y_dim-self.molecule.molsize2)/2)+")'>'"
			colors = {"P":"#B36AE2","T":"green"}
			for contact in self.pi_contacts.pi_contacts_for_drawing:
				ligand_atom_coords = []
				residue = str(contact[3]),str(contact[2]),str(contact[4])
				for atom_id in contact[0][0]:
					name = [atom.name  for atom in self.topology_data.universe.ligand.atoms if atom.id == atom_id][0]
					coord = self.molecule.ligand_atom_coords_from_diagr[name]
					ligand_atom_coords.append(coord)
				a = geometry.MultiPoint(ligand_atom_coords).centroid.coords.xy
				self.draw_pi_lines=self.draw_pi_lines+"<line x1='"+str(int(self.molecule.nearest_points_coords[residue][0]))+"' y1='"+str(int(self.molecule.nearest_points_coords[residue][1]))+"' x2='"+str(float(a[0][0]))+"' y2='"+str(float(a[1][0]))+"' style='stroke:"+colors[contact[1]]+";stroke-width:5' />"
				self.draw_pi_lines = self.draw_pi_lines+"<circle cx='"+str(int(a[0][0]))+"' cy='"+str(int(a[1][0]))+"' r='15' fill='"+colors[contact[1]]+"' />"
			self.draw_pi_lines = self.draw_pi_lines+"</g>"

	def draw_salt_bridges(self,color="blue"):
		"""
		For each bond that has been determined to be important, a line gets drawn.
		"""
		self.draw_saltbridges ="<g transform='translate("+str((self.molecule.x_dim-self.molecule.molsize1)/2)+","+str((self.molecule.y_dim-self.molecule.molsize2)/2)+")'>'"

		if self.saltbridges!=None:
			for bond in self.saltbridges.saltbridges_for_drawing:
				atom = self.topology_data.universe.atoms[bond[0]-1] #zero-based index vs one-based index
				residue = (atom.resname, str(atom.resid), atom.segid)
				self.draw_saltbridges=self.draw_saltbriges+"<line x1='"+str(int(self.molecule.nearest_points_coords[residue][0]))+"' y1='"+str(int(self.molecule.nearest_points_coords[residue][1]))+"' x2='"+str(float(self.molecule.ligand_atom_coords_from_diagr[bond[1]][0]))+"' y2='"+str(float(self.molecule.ligand_atom_coords_from_diagr[bond[1]][1]))+"' style='stroke:"+color+";stroke-width:4' />"
		self.draw_saltbridges= self.draw_saltbridges+"</g>"

	def add_smiles_id(self):
		for atom in self.lig_descr.ligand_atoms:
			self.lig_descr.ligand_atoms[atom]["SMILES ID"] = [v for k,v in self.molecule.atom_identities.items()if int(k)==atom][0]

	def normalise_colour(self,data_type="logP",color_type="viridis"):
		cmap = plt.get_cmap(color_type)
		min_value = min([v[data_type] for k,v in self.lig_descr.ligand_atoms.items()])
		max_value = max([v[data_type] for k,v in self.lig_descr.ligand_atoms.items()])
		norm = matplotlib.colors.Normalize(vmin=float("{0:.1f}".format(min_value)), vmax=float("{0:.1f}".format(max_value)))
		if data_type!="SASA" :
			color = {v["SMILES ID"]:cmap(norm(v[data_type])) for k,v in self.lig_descr.ligand_atoms.items()}
		else:
			color = {k:cmap(norm(v[data_type])) for k,v in self.lig_descr.ligand_atoms.items()}

		return color

	def normalise_size(self, data_type="MR"):
		min_value = min([v[data_type] for k,v in self.lig_descr.ligand_atoms.items()])
		max_value = max([v[data_type]+abs(min_value) for k,v in self.lig_descr.ligand_atoms.items()])
		if data_type!="SASA" :
			size = {v["SMILES ID"]:((v[data_type]+abs(min_value))/max_value)*0.4+0.2 for k,v in self.lig_descr.ligand_atoms.items()}
		else:
			size = {k:((v[data_type]+abs(min_value))/max_value)*0.4+0.2 for k,v in self.lig_descr.ligand_atoms.items()}
		return size

	def make_clouds(self,buff=90):
		#get fraction of size 
		buff = int(90 * float(self.molecule.molsize1)/900)
		polygons = [geometry.Point(point).buffer(buff) for point in self.molecule.ligand_atom_coords_from_diagr.values()]
		a =cascaded_union(polygons)

		self.shared_coords_x={}
		self.shared_coords_y={}

		for atom in self.molecule.ligand_atom_coords_from_diagr:
			point_coords=geometry.Point(self.molecule.ligand_atom_coords_from_diagr[atom]).buffer(buff).boundary.coords
			#get x atoms
			x_coords = [x for x in point_coords.xy[0]]
			y_coords = [x for x in point_coords.xy[1]]
			self.shared_coords_x[atom]=[x for x in a.exterior.coords.xy[0] if x in x_coords]
			self.shared_coords_y[atom]=[x for x in a.exterior.coords.xy[1] if x in y_coords]

	def draw_clouds(self,data_type="SASA",color_type="Blues"):

		fig_size_in_points_x = self.molecule.molsize1+int(200 * float(self.molecule.molsize1)/900) #To give enough expansion to the cloud - might have to update later
		fig_size_in_points_y = self.molecule.molsize2+int(200 * float(self.molecule.molsize1)/900) #To give enough expansion to the cloud - might have to update later
		color = self.normalise_colour(data_type,color_type)
		plt.figure(figsize=(float(fig_size_in_points_x)/72,float(fig_size_in_points_y)/72))
		
		#find if coordinates are under 0 and therefore will not show
		min_x = min([x for xds in self.shared_coords_x.values() for x in xds])
		if min_x >0:
			min_x = 0
		min_y = min([y for yds in self.shared_coords_y.values() for y in yds])
		if min_y >0:
			min_y =0 
		for atom in self.molecule.ligand_atom_coords_from_diagr:
			xs = self.shared_coords_x[atom]
			ys= self.shared_coords_y[atom]
			if len(xs) == 0:
				continue

			segs = [([xs[0]-min_x],[ys[0]-min_y])]
			for idx1, idx2 in zip(xrange(0, len(xs)-1), xrange(1, len(xs))):
				dx = xs[idx1] - xs[idx2]
				dy = ys[idx1] - ys[idx2]
				dist = (dx ** 2 + dy ** 2) ** 0.5
				if dist > 40:
					segs.append(([],[]))
				segs[-1][0].append(xs[idx2]-min_x)
				segs[-1][1].append(ys[idx2]-min_y)
			for seg_x, seg_y in segs:
				plt.plot(seg_x,seg_y,linewidth=8,c=color[[k for k,v in self.lig_descr.ligand_atoms.items() if v["name"]==atom][0]],transform=IdentityTransform()) 
		plt.axis('equal')
		plt.axis("off")


		pylab.savefig("cloud.svg",dpi=(100),transparent=True)
		self.manage_cloud_diagrams(min_x,min_y)

	def manage_cloud_diagrams(self,min_x,min_y):
		stop_line=None
		for i, line in enumerate(fileinput.input("cloud.svg", inplace=1)):
				if i <= 9 :
					continue
				if i > 10 and "<defs>" in line:
					stop_line=i
					continue
				if i>stop_line and stop_line!=None:
					continue
				else:
					sys.stdout.write(line)
		input1 =  '<g id="figure_1"'
		y_dim_cloud = self.molecule.molsize2+int(200 * float(self.molecule.molsize1)/900)
		output1 = "<g id='figure_1' transform='translate("+str((self.molecule.x_dim-self.molecule.molsize1)/2+min_x)+","+str((self.molecule.y_dim-self.molecule.molsize2)/2+y_dim_cloud+min_y)+") scale(1,-1)'>'"

		self.change_lines_in_svg('cloud.svg', input1, output1)
		input2 = "stroke-linecap:square;stroke-width:5;"
		output2 = "stroke-linecap:round;stroke-width:5;"
		self.change_lines_in_svg('cloud.svg', input2, output2)
		with open("cloud.svg", "r") as f:
			lines = f.readlines()
			self.cloud = self.cloud +"".join(map(str,lines))
			f.close()
	def get_rot_bonds(self):
		bond_ids=[]
		for atoms in self.lig_descr.rot_bonds:
			bond = self.topology_data.mol2.GetBondBetweenAtoms(atoms[0],atoms[1])
			bond_ids.append(bond.GetIdx())
		return bond_ids

	def draw_molecule(self,data_for_color,data_for_size,data_for_clouds,rot_bonds,color_for_clouds="Blues",color_type_color="viridis"):
		if data_for_color!=None:
			atom_colors = self.normalise_colour(data_for_color,color_type_color)
		else: 
			atom_colors = {}
		if data_for_size!=None:
			atom_size = self.normalise_size(data_for_size)
			if data_for_color==None:
				atom_colors = {x:(0.5,0.5,0.5) for x in range(self.topology_data.universe.ligand_noH.n_atoms)}
		else:
			atom_size = {}
		if rot_bonds!=None:
			bonds = self.get_rot_bonds()
		else:
			bonds = []
		self.molecule.load_molecule_in_rdkit_smiles((self.molecule.molsize1,self.molecule.molsize2),kekulize=True,bonds=bonds,bond_color=None,atom_color = atom_colors, size= atom_size )

		if data_for_clouds!=None:
			self.make_clouds()
			self.draw_clouds(data_for_clouds,color_for_clouds)
	def change_lines_in_svg(self,filename, string1,string2):
		"""Used to change lines in an SVG file. String1 is input - what is already present
		in the file, while string2 is output - the desired line. The input has to be searchable
		since .replace function is used. 
		Takes:
			* filename * - name of the SVG file to change
			* string1 * - input string
			* string2 * - output string

		"""
		for i,line in enumerate(fileinput.input(filename, inplace=1)):
			sys.stdout.write(line.replace(str(string1),str(string2)))


