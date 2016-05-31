__author__ = 'tparamo'

import numpy
from operator import itemgetter
import copy

valence = {"H":1, "C":4, "O":2, "N": 3, "S": 2, "P":5, "F":1,"B":3, "I":1, "BR":1, "SI":4, "CL": 1}
#Bond distances from Pyykko and Atsumi 2009
single_bond = {"H": 0.38, "C": 0.75, "O": 0.63, "N": 0.71, "S": 1.03, "P": 1.11, "F": 0.64 , "B": 0.85, "I": 1.33, "BR": 1.14, "SI": 1.16, "CL": 0.99}
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
    element = name[0]
    if len(name) > 1:
        key = str(name[:2]).upper()
        if (key in valence):
            element = str(name[:2]).upper()
    return element


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

                    # For aromatic rings
                    if get_atom_name(coordinates[i].name) == "C" and get_atom_name(coordinates[j].name) == "C" and not (coordinates[i].aromatic or coordinates[j].aromatic) and dist[i][j] < 1.40:
                        double = True
                        for k in range(j, len(coordinates)-1):
                            if k!=i and get_atom_name(coordinates[k].name) != "H" and dist[i][k]<dist[i][j]:
                                double = False
                        if double:
                            coordinates[i].num_bonds = coordinates[i].num_bonds + 1
                            coordinates[j].num_bonds = coordinates[j].num_bonds + 1
                            coordinates[i].aromatic = True
                            coordinates[j].aromatic = True
                            bond.type = 2
                    elif dist[i][j]>= (db - 0.02) and  dist[i][j]<= (db + 0.02) and coordinates[i].num_bonds<=(valence[get_atom_name(coordinates[i].name)]-1) and coordinates[j].num_bonds<=(valence[get_atom_name(coordinates[j].name)]-1):
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

                    coordinates[i].bonds.append(bond)

                    #Debug
                    #print coordinates[i].name + " " + coordinates[j].name + " " + str(dist[i][j]) + " " + str(db) + " " + str(coordinates[i].num_bonds) + " " + str(coordinates[j].num_bonds) + " " + str(bond.type)

    return num_bonds


def get_all_bonds(coordinates, name):
    bonds = []
    for i in range(0, len(coordinates)):
        for bond in coordinates[i].bonds:
            if coordinates[bond.a].name == name or coordinates[bond.b].name == name:
                bonds = bonds + [bond]
    return bonds



def reassign_bond(coordinates, bond, type):
    correction = type - bond.type
    bond.type = type
    coordinates[bond.a].num_bonds = coordinates[bond.a].num_bonds + correction
    coordinates[bond.a].charge = coordinates[bond.a].num_bonds - valence[get_atom_name(coordinates[bond.a].name)]
    coordinates[bond.b].num_bonds = coordinates[bond.b].num_bonds + correction
    coordinates[bond.b].charge = coordinates[bond.b].num_bonds - valence[get_atom_name(coordinates[bond.b].name)]


def sanitise_charge(coordinates, dist):
    #Final sanitisation - add more rules id there are more you can think of...

    correct = True

    for i in range(0, len(coordinates)):

        if get_atom_name(coordinates[i].name) == "P" and coordinates[i].charge == 1:
            #phosphate
            correct = False
            for bond in get_all_bonds(coordinates, coordinates[i].name):
                if bond.type==2:
                    reassign_bond(coordinates, bond, 1)
                    correct = True
                    break

        elif get_atom_name(coordinates[i].name) == "S" and coordinates[i].charge>0:
            #sulphate: I migh need other one or two double bonds
            correct = False
            cont = 2

            dists = []
            for bond in get_all_bonds(coordinates, coordinates[i].name):
                if (get_atom_name(coordinates[bond.a].name) == "O" and coordinates[bond.a].charge == -1) or (get_atom_name(coordinates[bond.b].name) == "O" and coordinates[bond.b].charge == -1):
                    dists.append((bond, dist[bond.a, bond.b]))

            if len(dists)>=1:
                sorted(dists, key=itemgetter(1))
                cont = 2
                for bond,d in dists:
                    if bond.type == 1 and cont>0:
                        reassign_bond(coordinates, bond, 2)
                        cont = cont - 1

                if cont==0:
                    #Valence is 6
                    coordinates[i].charge = coordinates[i].num_bonds - 6
                    correct = True
                elif cont==1:
                    #Valence is 4
                    coordinates[i].charge = coordinates[i].num_bonds - 4
                    correct = True

        elif get_atom_name(coordinates[i].name) == "N" and coordinates[i].charge == 0:
            # nitrate
            charge = 0
            bonds = []
            for bond in get_all_bonds(coordinates, coordinates[i].name):
                if get_atom_name(coordinates[bond.a].name) == "O":
                    charge = charge + coordinates[bond.a].charge
                    bonds = bonds + [bond]
                elif get_atom_name(coordinates[bond.b].name) == "O":
                    charge = charge + coordinates[bond.b].charge
                    bonds = bonds + [bond]

            if charge==-2:
                correct = False
                for bond in bonds:
                    reassign_bond(coordinates, bond, 2)
                    correct = True
                    break

        elif get_atom_name(coordinates[i].name) == "C" and coordinates[i].charge == 1:
            correct = False
            dists = []
            for bond in get_all_bonds(coordinates, coordinates[i].name):
                if (get_atom_name(coordinates[bond.a].name) in ["N", "C"] and coordinates[bond.a].charge == 1) or (get_atom_name(coordinates[bond.b].name) in ["N", "C"] and coordinates[bond.b].charge == 1):
                    dists.append((bond, dist[bond.a, bond.b]))

            if len(dists) >= 1:
                sorted(dists, key=itemgetter(1), reverse=True)
                reassign_bond(coordinates, dists[0][0], 1)
                correct = True

        elif get_atom_name(coordinates[i].name)=="C" and coordinates[i].charge == -1:
            dists = []
            aromatic =  False

            for bond in get_all_bonds(coordinates, coordinates[i].name):
                if coordinates[bond.a].aromatic or coordinates[bond.b].aromatic:
                    aromatic = True
                    break
            if not aromatic:
                for bond in get_all_bonds(coordinates, coordinates[i].name):
                    if (get_atom_name(coordinates[bond.a].name) == "N" and coordinates[bond.a].charge <= 0) or (get_atom_name(coordinates[bond.b].name) == "N" and coordinates[bond.b].charge <= 0):
                        dists.append((bond, dist[bond.a, bond.b]))

                if len(dists) >= 1:
                    sorted(dists, key=itemgetter(1))
                    reassign_bond(coordinates, dists[0][0], 2)


        if correct==False:
            raise AssertionError("Could not determine the structure correctly.")


def reset_ring(coordinates):
    ring_bonds = []

    for i in range(0, len(coordinates)):
        for bond in coordinates[i].bonds:
            if get_atom_name(coordinates[bond.a].name) in ["C","N"] and get_atom_name(coordinates[bond.b].name) in ["C","N"]:
                if (coordinates[bond.a].aromatic or coordinates[bond.a].charge != 0) or (coordinates[bond.b].aromatic or coordinates[bond.b].charge != 0):
                    ring_bonds = ring_bonds + [bond]


    # purge non-cyclic bonds
    for bond in ring_bonds:
        cont_i = 0
        cont_j = 0
        for subbond in ring_bonds:
            if coordinates[bond.a].name in [coordinates[subbond.a].name, coordinates[subbond.b].name]:
                cont_i = cont_i + 1
            if coordinates[bond.b].name in [coordinates[subbond.a].name, coordinates[subbond.b].name]:
                cont_j = cont_j + 1
        if not (cont_i>=2 and cont_j>=2):
            ring_bonds.remove(bond)

    for bond in ring_bonds:
        neg_charge = coordinates[bond.a].charge == -1 or coordinates[bond.b].charge == -1
        #pos_charge = (get_atom_name(coordinates[bond.a].name) == "C" and coordinates[bond.a].charge == 1) or (get_atom_name(coordinates[bond.b].name) == "C" and coordinates[bond.b].charge == 1)
        if neg_charge: #or pos_charge:

            start = True
            stop = False
            current_bond = bond
            odd = True

            ring_bonds_buffer = ring_bonds

            while (current_bond.id != bond.id or start) and not stop:
                #print coordinates[current_bond.a].name + " " + coordinates[current_bond.b].name
                start = False
                if odd:
                    reassign_bond(coordinates, current_bond, 2)
                else:
                    reassign_bond(coordinates, current_bond, 1)

                odd = not odd

                stop = True
                candidate = ""

                for next_bond in ring_bonds_buffer:
                    if current_bond.id != next_bond.id:
                        righthand = (coordinates[current_bond.b].name == coordinates[next_bond.a].name)
                        lefthand = (coordinates[current_bond.b].name == coordinates[next_bond.b].name)
                        if righthand:
                            current_bond = next_bond
                            stop = False
                            break
                        elif lefthand:
                            candidate = next_bond

                if stop==True and candidate!="":
                    a = candidate.a
                    b = candidate.b
                    current_bond = candidate
                    current_bond.a = b
                    current_bond.b = a
                    stop = False

                if not stop:
                    ring_bonds_buffer.remove(current_bond)


def sanitise_aromaticity(coordinates, dist):

    for i in range(0, len(coordinates)):
        if get_atom_name(coordinates[i].name) in ["C","N"] and coordinates[i].charge == -1:
            min_dist = 100000
            min_bond = []

            for bond in get_all_bonds(coordinates, coordinates[i].name):
                lefthand = get_atom_name(coordinates[bond.a].name) in ["C","N","O"] and coordinates[bond.a].charge == -1
                righthand = get_atom_name(coordinates[bond.b].name) in ["C","N","O"] and coordinates[bond.b].charge == -1
                if lefthand and righthand:
                    if dist[bond.a, bond.b] < min_dist:
                        min_bond = bond
                        min_dist = dist[bond.a, bond.b]

            if min_bond != []:
                reassign_bond(coordinates, min_bond, 2)



def calculate_formal_charge(coordinates, dist):
    for i in range(0, len(coordinates)):
        if coordinates[i].num_bonds!=valence[get_atom_name(coordinates[i].name)]:
            coordinates[i].charge = coordinates[i].num_bonds - valence[get_atom_name(coordinates[i].name)]
    #Sanitise non standard bonds
    sanitise_aromaticity(coordinates, dist)
    sanitise_charge(coordinates, dist)

    for i in range(0, len(coordinates)):
        if get_atom_name(coordinates[i].name) in ["C", "N"] and coordinates[i].charge == -1:
            reset_ring(coordinates)


def read_PDB(pdb):
    coordinates = []
    try:
        pdb_file = open(pdb, "r")
        for line in pdb_file.readlines():
            if line.startswith("ATOM") or line.startswith("HETATM"):
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
            element = get_atom_name(coordinates[i].name)
            if len(element)==2:
                element = element[0]+element[1].lower()

            atom_info = ' {:6d}'.format(i + 1) + '  {:7}'.format(coordinates[i].name) + ' {:9.4f}'.format(
                coordinates[i].x) + ' {:9.4f}'.format(coordinates[i].y) + ' {:9.4f}'.format(coordinates[i].z)
            atom_info = atom_info + ' {:5}'.format(element) + ' {:4}'.format(
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

if __name__ == '__main__':
    import sys
    pdb2mol2(sys.argv[1])



