__author__ = 'tparamo'

import numpy


valence = {"H":1, "C":4, "O":2, "N": 3, "S": 2, "P":5, "F":1,"B":3, "I":1, "BR":1, "SI":4, "CL": 1}
#Bond distances from Pyykko and Atsumi 2009
single_bond = {"H": 0.32, "C": 0.75, "O": 0.63, "N": 0.71, "S": 1.03, "P": 1.11, "F": 0.64 , "B": 0.85, "I": 1.33, "BR": 1.14, "SI": 1.16, "CL": 0.99}
double_bond = {"H":0.32, "C":0.67, "O":0.57, "N": 0.60, "S": 0.94, "P":1.02, "F":0.59,"B":0.78, "I":1.29, "BR":1.09, "SI":1.07, "CL": 0.95}


class Atom(object):
    def __init__(self, name, x, y, z, resid, resname):
        self.name = name
        self.x = x
        self.y = y
        self.z = z
        self.resid = resid
        self.resname = resname
        self.num_bonds = 0
        self.bonds = []
        self.charge = 0
        self.aromatic = False


class Bond(object):
    def __init__(self, id, a, b, type):
        self.id = id
        self.a = a
        self.b = b
        self.type = type


def get_atom_name(name):
    name = name[0]
    if len(name) > 2:
        key = str(name[:2]).upper()
        if (key in vdw):
            name = name[name[:2]]
    return name


def distance(atom_a, atom_b):
    sqr = ((atom_a.x - atom_b.x) ** 2.0) + ((atom_a.y - atom_b.y) ** 2.0) + ((atom_a.z - atom_b.z) ** 2.0)
    return sqr ** 0.5


def calculate_distance_matrix(coordinates):
    num_residues = len(coordinates)
    matrix = numpy.zeros((num_residues, num_residues))
    for i in range(0, num_residues):
        for j in range(0, num_residues):
            matrix[i][j] = distance(coordinates[i], coordinates[j])
    return matrix


def calculate_bonds(coordinates, dist):
    num_bonds = 0
    for i in range(0, len(coordinates)):
        for j in range(0,len(coordinates)):
            if i<j:
                sb = single_bond[get_atom_name(coordinates[i].name)] + single_bond[get_atom_name(coordinates[j].name)]
                #Debug
                #print coordinates[i].name + " " + coordinates[j].name + " " + str(dist[i][j]) + " " + str(sb)
                if dist[i][j]<(sb + 0.1):
                    num_bonds = num_bonds + 1
                    coordinates[i].num_bonds = coordinates[i].num_bonds + 1
                    coordinates[j].num_bonds = coordinates[j].num_bonds + 1
                    bond = Bond(num_bonds, i,j, 1)

                    #Determine the type of bond based on distance
                    db = double_bond[get_atom_name(coordinates[i].name)] + double_bond[get_atom_name(coordinates[j].name)]

                    if dist[i][j]>= (db - 0.02) and  dist[i][j]<= (db + 0.02) and coordinates[i].num_bonds<=(valence[get_atom_name(coordinates[i].name)]-1) and coordinates[j].num_bonds<=(valence[get_atom_name(coordinates[j].name)]-1):
                        coordinates[i].num_bonds = coordinates[i].num_bonds + 1
                        coordinates[j].num_bonds = coordinates[j].num_bonds + 1
                        bond.type = 2
                    elif dist[i][j]< (db - 0.02):
                        if coordinates[i].num_bonds<=(valence[get_atom_name(coordinates[i].name)]-2) and coordinates[j].num_bonds<=(valence[get_atom_name(coordinates[j].name)]-2):
                            coordinates[i].num_bonds = coordinates[i].num_bonds + 2
                            coordinates[j].num_bonds = coordinates[j].num_bonds + 2
                            bond.type = 3
                        elif coordinates[i].num_bonds<=(valence[get_atom_name(coordinates[i].name)]-1) and coordinates[j].num_bonds<=(valence[get_atom_name(coordinates[j].name)]-1):
                            coordinates[i].num_bonds = coordinates[i].num_bonds + 1
                            coordinates[j].num_bonds = coordinates[j].num_bonds + 1
                            bond.type = 2
                    #For aromatic rings
                    elif get_atom_name(coordinates[i].name)=="C" and get_atom_name(coordinates[j].name)=="C" and not (coordinates[i].aromatic or coordinates[j].aromatic) and dist[i][j]<1.40:
                        coordinates[i].num_bonds = coordinates[i].num_bonds + 1
                        coordinates[j].num_bonds = coordinates[j].num_bonds + 1
                        coordinates[i].aromatic = True
                        coordinates[j].aromatic = True
                        bond.type = 2

                    coordinates[i].bonds.append(bond)

                    #Debug
                    #print coordinates[i].name + " " + coordinates[j].name + " " + str(dist[i][j]) + " " + str(db) + " " + str(coordinates[i].num_bonds) + " " + str(coordinates[j].num_bonds)

    return num_bonds


def sanitise_charge(coordinates, dist):
    #Final sanitisation - add more rules id there is more you can think of...
    for i in range(0, len(coordinates)):
        if get_atom_name(coordinates[i].name) == "C" and coordinates[i].charge == -1:
            min_dist = 100000
            min_bond = -1
            for index,bond in enumerate(coordinates[i].bonds):
                if get_atom_name(coordinates[bond.b].name) == "N" and coordinates[bond.b].charge==0 and coordinates[i].charge == -1:
                    # Is available to be charged
                    if dist[bond.a, bond.b]<min_dist:
                        min_bond = index
                        min_dist = dist[bond.a, bond.b]

                if min_bond>=0:
                    coordinates[i].bonds[min_bond].type = 2
                    coordinates[i].num_bonds = coordinates[i].num_bonds + 1
                    coordinates[i].charge = coordinates[i].num_bonds - valence[get_atom_name(coordinates[i].name)]
                    coordinates[coordinates[i].bonds[min_bond].b].num_bonds = coordinates[coordinates[i].bonds[min_bond].b].num_bonds + 1
                    coordinates[coordinates[i].bonds[min_bond].b].charge = coordinates[coordinates[i].bonds[min_bond].b].num_bonds - valence[get_atom_name(coordinates[coordinates[i].bonds[min_bond].b].name)]
                else:
                    print "Could not determine the structure correctly- check your input file?\n"


def sanitise_aromatic_ring(coordinates):
    for i in range(0, len(coordinates)):
        if get_atom_name(coordinates[i].name) == "C" and coordinates[i].charge == -1:
            for bond in coordinates[i].bonds:
                candidate = bond.b
                if get_atom_name(coordinates[candidate].name) == "C" and coordinates[candidate].charge == -1 and coordinates[i].charge == -1:
                    bond.type = 2
                    coordinates[i].num_bonds = coordinates[i].num_bonds + 1
                    coordinates[i].charge = coordinates[i].num_bonds - valence[get_atom_name(coordinates[i].name)]
                    coordinates[candidate].num_bonds = coordinates[candidate].num_bonds + 1
                    coordinates[candidate].charge = coordinates[candidate].num_bonds - valence[get_atom_name(coordinates[candidate].name)]
                    break


def calculate_formal_charge(coordinates, dist):
    for i in range(0, len(coordinates)):
        if coordinates[i].num_bonds!=valence[get_atom_name(coordinates[i].name)]:
            coordinates[i].charge = coordinates[i].num_bonds - valence[get_atom_name(coordinates[i].name)]
    #Sanitise poorly defined aromatic rings and charges
    sanitise_aromatic_ring(coordinates)
    sanitise_charge(coordinates, dist)


def read_PDB(pdb):
    coordinates = []
    try:
        pdb_file = open(pdb, "r")
        for line in pdb_file.readlines():
            if line.__contains__("ATOM") or line.__contains__("HETATM"):
                atom = str(line[12:16]).replace(" ","")
                x = float(line[31:38])
                y = float(line[39:46])
                z = float(line[47:54])
                resid = line[23:26]
                resname = line[17:20]

                coordinates.append(Atom(atom, x, y, z, resid, resname))

        pdb_file.close()

    except IOError:
        print pdb + "does not exist\n"

    return coordinates


def write_mol2(pdb, num_bonds, coordinates):
    try:
        mol2 = open(pdb[:-4]+"_test.mol2", "w")
        mol2.writelines(["@<TRIPOS>MOLECULE\n", pdb+"\n", " %2i %2i %2i %2i %2i\n" % (len(coordinates), num_bonds, 0,0,0), "SMALL\n", "FORMAL\n", "\n", "@<TRIPOS>ATOM\n"])

        for i in range(0, len(coordinates)):
            atom_info = ' {:6d}'.format(i + 1) + '  {:7}'.format(coordinates[i].name) + ' {:9.4f}'.format(
                coordinates[i].x) + ' {:9.4f}'.format(coordinates[i].y) + ' {:9.4f}'.format(coordinates[i].z)
            atom_info = atom_info + ' {:5}'.format(get_atom_name(coordinates[i].name)) + ' {:4}'.format(
                coordinates[i].resid) + ' {:6}'.format(coordinates[i].resname + coordinates[i].resid) + ' {:6d}'.format(
                coordinates[i].charge)
            mol2.write(atom_info + "\n")

        mol2.write("@<TRIPOS>BOND\n")

        for i in range(0, len(coordinates)):
            for bond in coordinates[i].bonds:
                mol2.write(" %5i %5i %5i %4i\n" % (bond.id, bond.a+1, bond.b+1, bond.type))
        mol2.close()

    except IOError:
        print("Can't write!\n")


def write_pdb(pdb, num_bonds, coordinates):
    try:
        pdb = open(pdb[:-4]+"_test.pdb", "w")

        for i in range(0, len(coordinates)):
            atom_info = "ATOM  " + "{:>5d}".format(i + 1) + "  {:3}".format(coordinates[i].name) + " {:4}".format(coordinates[i].resname) + " {:4}".format(coordinates[i].resid)
            atom_info = atom_info + "    {:8.3f}".format(coordinates[i].x) + "{:8.3f}".format(coordinates[i].y) + "{:8.3f}".format(coordinates[i].z) + "  1.00  0.00"
            atom_info = atom_info + "{:>12}".format(get_atom_name(coordinates[i].name))
            pdb.write(atom_info + "\n")

        for i in range(0, len(coordinates)):
            for bond in coordinates[i].bonds:
                pdb.write("CONECT %4i %4i %4i\n" % (bond.id, bond.a+1, bond.b+1))

        pdb.write("END\n")
        pdb.close()

    except IOError:
        print("Can't write!\n")

def pdb2mol2(pdb, mol2=True):
    coordinates = read_PDB(pdb)
    if len(coordinates) > 0:
        dist = calculate_distance_matrix(coordinates)
        num_bonds = calculate_bonds(coordinates, dist)
        calculate_formal_charge(coordinates, dist)
        if mol2:
            write_mol2(pdb, num_bonds, coordinates)
        else:
            write_pdb(pdb, num_bonds, coordinates)
    else:
        print "There are no atoms!\n"


