import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from IPython.display import SVG
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdDepictor
from shapely import geometry
import numpy as np
from itertools import combinations
import colorsys

class Molecule(object):
    def __init__(self, mol2_file, topol_object, rmsf_object=None):
        self.svg = None
        self.universe = topol_object
        self.rmsf = rmsf_object
        self.final_svg = None
        self.load_molecule_in_rdkit(mol2_file)
        self.atom_coords_from_diagramm = {}
        self.nearest_points ={}
        self.nearest_points_projection = {}
        self.nearest_points_coords ={}
        self.coefficient ={}
        self.convex_hull()
        self.make_new_projection_values()
    def pseudocolor(self,val, minval, maxval):
        # convert val in range minval..maxval to the range 0..120 degrees which
        # correspond to the colors red..green in the HSV colorspace
        h = (float(val-minval) / (maxval-minval)) * 120
        # convert hsv color (h,1,1) to its rgb equivalent
        # note: the hsv_to_rgb() function expects h to be in the range 0..1 not 0..360
        r, g, b = colorsys.hsv_to_rgb(h/360, 1., 1.)
        return (r, g, b)
    def load_molecule_in_rdkit(self, mol2_file, molSize=(600,300),kekulize=True):
        self.ligand_in_rdkit=Chem.MolFromMol2File(mol2_file)
        highlight = []
        colors = {}
        rdDepictor.Compute2DCoords(self.ligand_in_rdkit)
        #rdDepictor.Compute2DCoords(self.ligand_in_rdkit)
        if self.rmsf is not None:
            for i in range(self.ligand_in_rdkit.GetNumAtoms()):
                highlight.append(i)
                colors[i] = self.pseudocolor(self.rmsf.ligand_rmsf[i], self.rmsf.min_value, self.rmsf.max_value)
        else:
            for i in range(self.ligand_in_rdkit.GetNumAtoms()):
                highlight.append(i)
                colors[i]=(1,1,1)
        drawer = rdMolDraw2D.MolDraw2DSVG(molSize[0],molSize[1])
        drawer.DrawMolecule(self.ligand_in_rdkit,highlightAtoms=highlight,highlightBonds=[], highlightAtomColors=colors)
        drawer.FinishDrawing()
        self.svg = drawer.GetDrawingText().replace('svg:','')
        filesvg = open("molecule.svg", "w+")
        filesvg.write(self.svg)
    
    def convex_hull(self):
        """Draws a convex hull around ligand atoms and expands it, giving space to put diagramms on"""
        #Get coordinates of ligand atoms (needed to draw the convex hull around)
        ligand_atom_coords = []
        with open ("molecule.svg", "r") as f:
            lines = f.readlines()
            i=0
            for line in lines:
                if line.startswith("<ellipse"): 
                    ligand_atom_coords.append([float(line.rsplit("'",10)[1]), float(line.rsplit("'",10)[3])]) 
                    for atom in range(len(self.universe.closest_atoms.keys())):
                        if self.universe.closest_atoms.values()[atom][1]==i:
                            self.atom_coords_from_diagramm[self.universe.closest_atoms.keys()[atom]] = [float(line.rsplit("'",10)[1]), float(line.rsplit("'",10)[3]),self.universe.closest_atoms.values()[atom][0]]
                    i+=1
        ligand_atom_coords=np.array(ligand_atom_coords)  
        # Get the convex hull around ligand atoms 
        self.a = geometry.MultiPoint(ligand_atom_coords).convex_hull

        self.b_for_all = {}
        for residue in self.atom_coords_from_diagramm:
            self.b_lenght = None
            b = self.a.boundary.parallel_offset(80,"left",join_style=2).convex_hull
            point =geometry.Point((self.atom_coords_from_diagramm[residue][0],self.atom_coords_from_diagramm[residue][1]))
            self.b_for_all[residue] = (b.boundary.project(point) % b.boundary.length) 
            self.b_lenght = b.boundary.length
        self.make_multiple_hulls()
        
        
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
            if coeff[index1] <= coeff[index2]:
                if self.b_lenght-coeff[index2]<self.b_lenght/5: #a quick and dirty solution, but works
                    forces[index1].append(f[1]) # push to left (smaller projection value) 
                    forces[index2].append(f[0])
                else:
                    #all is normal
                    forces[index1].append(f[0]) # push to left (smaller projection value) 
                    forces[index2].append(f[1])
            else:
                if self.b_lenght-coeff[index1]<self.b_lenght/5: #a quick and dirty solution, but works
                    forces[index1].append(f[0]) # push to left (smaller projection value) 
                    forces[index2].append(f[1])
                else:
                #if all is normal
                    forces[index1].append(f[1]) # push to left (smaller projection value) 
                    forces[index2].append(f[0])
        forces = {k:sum(v) for k,v in forces.items()}
        
        energy = sum([abs(x) for x in forces.values()])

        return [(forces[k]/10+v) for k, v in enumerate(values)], energy
  
    def make_new_projection_values(self):
        """Run do_step function until the diagramms have diverged from each other"""
        values = [v for v in self.nearest_points_projection.values()]
        xy_values = [v for v in self.nearest_points_coords.values()]
        coeff_value = [v for v in self.b_for_all.values()]
        #values = [x for x in startvalues]
        #xy_values = [x for x in startvalues1]
        energy = 100
        while energy > 0.2:
            values, energy = self.do_step(values,xy_values,coeff_value, width=90)
            i=0
            xy_values =[]
            for residue in  self.nearest_points_coords:
                b = self.a.boundary.parallel_offset(self.universe.closest_atoms[residue][2]*32.0+24,"left",join_style=2).convex_hull
                self.nearest_points_projection[residue] = values[i]
                self.nearest_points[residue] = b.boundary.interpolate(self.nearest_points_projection[residue] % b.boundary.length)
                self.nearest_points_coords[residue] = self.nearest_points[residue].x, self.nearest_points[residue].y
                xy_values.append(self.nearest_points_coords[residue])
                i+=1
            values = [v for v in self.nearest_points_projection.values()]
        self.x_dim  = max(x[0] for i,x in enumerate(xy_values))-min(x[0] for i,x in enumerate(xy_values))+216.00
        self.y_dim = max(x[1] for i,x in enumerate(xy_values))-min(x[1] for i,x in enumerate(xy_values))+216.00
        if self.x_dim<600:
            self.x_dim=600+216
        if self.y_dim<300:
            self.y_dim=300+216

    def make_multiple_hulls(self):
        for residue in self.atom_coords_from_diagramm:
            b = self.a.boundary.parallel_offset(self.universe.closest_atoms[residue][2]*32.0+24,"left",join_style=2).convex_hull
            point =geometry.Point((self.atom_coords_from_diagramm[residue][0],self.atom_coords_from_diagramm[residue][1]))
            self.nearest_points_projection[residue] = (b.boundary.project(point) % b.boundary.length)
            self.nearest_points[residue] = b.boundary.interpolate(self.nearest_points_projection[residue] % b.boundary.length)
            self.nearest_points_coords[residue]=self.nearest_points[residue].x,self.nearest_points[residue].y
    