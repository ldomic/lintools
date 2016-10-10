import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdDepictor
from shapely import geometry
import numpy as np
from itertools import combinations
import operator

class Molecule(object):
    """
    This class deals with the drawing of the 2D molecule in RDkit environment. Once the data has been imported and 2D
    representation drawn, the coordinates of each atom are obtained from the drawing. Since it is known, which atom 
    is close to each plotted residue, the residues are placed in vicinity of that atom and the overlap treated by 
    minimising distances in 2D. The resulting 2D coordinates where each residue should be placed in the final image 
    are inherited by Figure class.

    Takes:
        * topology_data_object * - information about the system (lintools.Data object)

    Initializing the object will lead to execution of the functions present in the class 
    providing necessary material for next steps, namely assembly of figure. This was done
    since it is very unlikely that the process would ever have to be done seperately step 
    by step. 
    """
    __version__= "09.2016"
    def __init__(self, topology_data_object):
        self.topology_data = topology_data_object
        self.ligand_atom_coords_from_diagr={}
        self.nearest_points ={}
        self.nearest_points_projection = {}
        self.nearest_points_coords ={}
        self.coefficient ={}
        self.ligand_atom_coords = []
        self.arc_coords=None
        self.load_molecule_in_rdkit_smiles()
        self.convex_hull()
        self.make_new_projection_values()
    def load_molecule_in_rdkit_smiles(self, molSize=(900,450),kekulize=True):
        """
        Loads mol2 file in rdkit without the hydrogens - they do not have to appear in the final
        figure. Once loaded, the molecule is converted to SMILES format which RDKit appears to 
        draw best - since we do not care about the actual coordinates of the original molecule, it
        is sufficient to have just 2D information. 

        Some molecules can be problematic to import and steps such as stopping sanitize function can
        be taken. This is done automatically if problems are observed. However, better solutions can
        also be implemented and need more research.

        The molecule is then drawn from SMILES in 2D representation without hydrogens. The drawing is 
        saved as an SVG file.
        """
        highlight=[]
        colors={}
        mol2_in_rdkit = self.topology_data.mol2 #need to reload without hydrogens
        try:
            mol2_in_rdkit = Chem.RemoveHs(mol2_in_rdkit)
            self.smiles = Chem.MolFromSmiles(Chem.MolToSmiles(mol2_in_rdkit))
        except ValueError:
            mol2_in_rdkit = Chem.RemoveHs(mol2_in_rdkit, sanitize = False)
            self.smiles = Chem.MolFromSmiles(Chem.MolToSmiles(mol2_in_rdkit), sanitize=False)
        self.atom_identities = {}
        i=0
        for atom in self.smiles.GetAtoms():
            self.atom_identities[mol2_in_rdkit.GetProp('_smilesAtomOutputOrder')[1:].rsplit(",")[i]] = atom.GetIdx()
            i+=1
        mc = Chem.Mol(self.smiles.ToBinary())
        if kekulize:
            try:
                Chem.Kekulize(mc)
            except:
                mc = Chem.Mol(self.smiles.ToBinary())
        if not mc.GetNumConformers():
            rdDepictor.Compute2DCoords(mc)
        for i in range(mol2_in_rdkit.GetNumAtoms()):
            highlight.append(i)
            colors[i]=(1,1,1)
        drawer = rdMolDraw2D.MolDraw2DSVG(molSize[0],molSize[1])
        drawer.DrawMolecule(mc,highlightAtoms=highlight,highlightBonds=[], highlightAtomColors=colors)
        drawer.FinishDrawing()
        self.svg = drawer.GetDrawingText().replace('svg:','')
        filesvg = open("molecule.svg", "w+")
        filesvg.write(self.svg)
    def convex_hull(self):
        """
        Draws a convex hull around ligand atoms and expands it, giving space to put diagramms on.
        This is done with the help of Shapely.geometry class. The initial convex hull the residue 
        coordinates are inserted on, determines the order the coordinates are going to be moved, i.e.
        if the residue 1 is on the right side of residue 2, it will be pushed to the right, while 
        residue 2 will be moved to the left.


        Also determines the 2D coordinates of all atoms in drawing and makes a list with those.

        """
        #Get coordinates of ligand atoms (needed to draw the convex hull around)
        
        ligand_atoms = [x.name for x in self.topology_data.universe.ligand_noH.atoms]
        with open ("molecule.svg", "r") as f:
            lines = f.readlines()
            i=0
            for line in lines:
                if line.startswith("<ellipse"): 
                    self.ligand_atom_coords.append([float(line.rsplit("'",10)[1]), float(line.rsplit("'",10)[3])]) 
                    for atom_id in self.atom_identities:
                        if i == self.atom_identities[atom_id]:
                            self.ligand_atom_coords_from_diagr[ligand_atoms[int(atom_id)]]=[float(line.rsplit("'",10)[1]), float(line.rsplit("'",10)[3])]
                    i+=1
                    
        self.ligand_atom_coords=np.array(self.ligand_atom_coords)  
        self.a = geometry.MultiPoint(self.ligand_atom_coords).convex_hull
        self.b = self.a.boundary.buffer(120).convex_hull
        self.b_for_all ={}
        self.b_lenght = self.b.boundary.length
        for residue in self.topology_data.closest_atoms:
            mean_distance =np.array([x[1] for x in self.topology_data.closest_atoms[residue]]).mean()
            b = self.a.boundary.parallel_offset(mean_distance*50+50,"left",join_style=2).convex_hull
            projection =[]
            projection_init = []
            for atom in self.topology_data.closest_atoms[residue]:
                point =geometry.Point((self.ligand_atom_coords_from_diagr[atom[0]][0],self.ligand_atom_coords_from_diagr[atom[0]][1]))
                projection.append(abs(b.boundary.project(point) % b.boundary.length))
                projection_init.append(abs(self.b.boundary.project(point) % self.b.boundary.length))
            # Check whether projections are not crossing the boundary point (i.e. end of circle) - is one number in the projection very different from any other?
            for (index1,number1), (index2,number2) in combinations(enumerate(projection),2):
                if abs(number1-number2)>b.boundary.length/2:
                    proj =[]
                    for atom in projection:
                        if atom == max([number1,number2]):
                            proj.append(atom-b.boundary.length)
                        else:
                            proj.append(atom)
                    projection = proj
            for (index1,number1), (index2,number2) in combinations(enumerate(projection_init),2):
                if abs(number1-number2)>self.b.boundary.length/2:
                    proj =[]
                    for atom in projection_init:
                        if atom == max([number1,number2]):
                            proj.append(atom-self.b.boundary.length)
                        else:
                            proj.append(atom)
                    projection_init = proj
            self.nearest_points_projection[residue] = np.array(projection).mean()
            self.b_for_all[residue] = np.array(projection_init).mean()
            self.nearest_points[residue] = b.boundary.interpolate(self.nearest_points_projection[residue] % b.boundary.length)
            self.nearest_points_coords[residue]=self.nearest_points[residue].x,self.nearest_points[residue].y


   

    def calc_2d_forces(self,x1,y1,x2,y2,width):
        """Calculate overlap in 2D space"""
        #calculate a
        if x1>x2:
            a = x1-x2
        else:
            a = x2-x1

        a_sq=a*a
        #calculate b
        if y1>y2:
            b = y1-y2
        else: 
            b = y2-y1

        b_sq=b*b
        
        #calculate c
        from math import sqrt
        c_sq = a_sq+b_sq

        c = sqrt(c_sq)

        if c > width:
            return 0,0
        else:
            overlap = width-c
        return -overlap/2, overlap/2

    
    def do_step(self, values, xy_values,coeff, width):
        """Calculates forces between two diagrams and pushes them apart by tenth of width"""
        forces = {k:[] for k,i in enumerate(xy_values)}
        for (index1, value1), (index2,value2) in combinations(enumerate(xy_values),2):
            f = self.calc_2d_forces(value1[0],value1[1],value2[0],value2[1],width)
            if coeff[index1] < coeff[index2]:
                if self.b_lenght-coeff[index2]<self.b_lenght/10: #a quick and dirty solution, but works
                    forces[index1].append(f[1]) # push to left (smaller projection value) 
                    forces[index2].append(f[0])
                else:
                    #all is normal
                    forces[index1].append(f[0]) # push to left (smaller projection value) 
                    forces[index2].append(f[1])
            else:
                if self.b_lenght-coeff[index1]<self.b_lenght/10: #a quick and dirty solution, but works
                    forces[index1].append(f[0]) # push to left (smaller projection value) 
                    forces[index2].append(f[1])
                else:
                #if all is normal
                    forces[index1].append(f[1]) # push to left (smaller projection value) 
                    forces[index2].append(f[0])
        forces = {k:sum(v) for k,v in forces.items()}
        
        energy = sum([abs(x) for x in forces.values()])

        return [(forces[k]/10+v) for k, v in enumerate(values)], energy
  
    

    def make_new_projection_values(self,width=160):
        """Run do_step function until the diagramms have diverged from each other.
        Also determines how big the figure is going to be by calculating the borders 
        from new residue coordinates. These are then added some buffer space.

        """
        #Make gap between residues bigger if plots have a lot of rings - each ring after the 4th
        #give extra 12.5px space
        if self.topology_data.ring_number>4:
            width = width + (self.topology_data.ring_number-4)*12.5
        values = [v for v in self.nearest_points_projection.values()]
        xy_values = [v for v in self.nearest_points_coords.values()]
        coeff_value = [v for v in self.b_for_all.values()]
        energy = 100
        while energy > 0.2:
            values, energy = self.do_step(values,xy_values,coeff_value, width)
            i=0
            xy_values =[]
            for residue in  self.nearest_points_coords:
                b = self.a.boundary.parallel_offset(self.topology_data.closest_atoms[residue][0][1]*50+50,"left",join_style=2).convex_hull
                self.nearest_points_projection[residue] = values[i]
                self.nearest_points[residue] = b.boundary.interpolate(self.nearest_points_projection[residue] % b.boundary.length)
                self.nearest_points_coords[residue] = self.nearest_points[residue].x, self.nearest_points[residue].y
                xy_values.append(self.nearest_points_coords[residue])
                i+=1
            values = [v for v in self.nearest_points_projection.values()]
        
        #Calculate the borders of the final image
        max_x = int(max(v[0] for k,v in self.nearest_points_coords.items()))
        min_x = int(min(v[0] for k,v in self.nearest_points_coords.items()))
        min_y = int(min(v[1] for k,v in self.nearest_points_coords.items()))
        max_y = int(max(v[1] for k,v in self.nearest_points_coords.items()))
        if min_x<0:
            self.x_dim =(max_x-min_x)+600 #600 acts as buffer
        elif max_x<900 and min_x<0: #In case all residues are grouped on one end of the molecule
            self.x_dim = (900-min_x)+600
        elif max_x<900 and min_x>0:
            self.x_dim = 900+600
        else:
            self.x_dim = max_x+600
        if min_y<0:
            self.y_dim = (max_y-min_y)+400 #400 acts as buffer
        elif max_y<450 and min_y<0:
            self.y_dim = (450-min_y)+400
        elif max_y<450 and min_y>0:
            self.y_dim = 450+400
        else:
            self.y_dim = max_y+400