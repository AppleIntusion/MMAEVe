''' ______________________________________________________________________
|                                                                    |
| Title:   Builder foundation                                        |
|                                                                    |
| Author:  Gubbin Eel (Samuel Lindsay, Independent Project)          |
|                                                                    |
| Purpose: Defines the Atom() class and all of the functions that    |
|          are used to build micelles, vesicles, and bilayers.       |
|                                                                    |
| Usage:   The other 'foundational' portions of the program will all |
|          contain the 'import foundation as fd' statement at their  |
|          start                                                     |
|____________________________________________________________________|
'''

'''--------
| Modules |
--------'''

import string
import copy
import numpy as np

'''-----------------
| Global Variables |
-----------------'''

####################################################################
# Atomic masses of the most abundant isotope of each of the most   #
# common elements in biological molecules                          #
#                                                                  #
# Mass data:                                                       #
#     G. Audi, A. H. Wapstra Nucl. Phys A. 1993, 565, 1-65         #
#     G. Audi, A. H. Wapstra Nucl. Phys A. 1995, 595, 409-480      #
#                                                                  #
# Abundance:                                                       #
#     K. J. R. Rosman, P. D. P. Taylor, Pure Appl. Chem. 1999, 71, #
#     1593-1607                                                    #
####################################################################
atomic_masses = {'H': 1.007825, 'He': 4.002603, 'Li': 7.016004, \
                 'Be': 9.012182, 'B': 11.009305, 'C': 12.000000, \
                 'N': 14.003074, 'O': 15.994915, 'F': 18.998403, \
                 'P': 30.973762, 'S': 31.972071}

'''--------
| Classes |
--------'''


class Atom(object):
    '''
    Represents an atom in 3D space

    Attributes: x, y, z, mass, name, resn
    '''

    def __init__(self, x=0, y=0, z=0, mass=0, name='', resn=''):
        self.x = x
        self.y = y
        self.z = z
        self.mass = mass
        self.name = name
        self.resn = resn

    def __sub__(self, other):
        x, y, z = self.x - other.x, self.y - other.y, self.z - other.z
        new = [x, y, z]
        return new

    # Purpose:   To translate the position of an 'Atom' object by
    #            the specified amount. A pure function; does not
    #            update the argument.
    # Arguments: (1) 'Atom' object; use method format
    #            (2) List or tuple, 3 elements long. The amount by
    #                which each coordinate will be translated
    # Modules:   copy
    # Functions: NONE
    def translate_atom(self, t=[0, 0, 0]):
        position = copy.deepcopy(self)
        position.x += t[0]
        position.y += t[1]
        position.z += t[2]
        return position

    # Purpose:   To calcualte the distance between two 'Atom'
    #            objects
    # Arguments: (1) 'Atom' object; use method format
    #            (2) 'Atom' object; provide as function argument
    # Modules:   NONE
    # Functions: NONE
    def distance(self, other):
        dist = (((self.x - other.x) ** 2) + \
                ((self.y - other.y) ** 2) + \
                ((self.z - other.z) ** 2)) ** (1 / 2)
        return dist

    # Purpose:   To rotate an 'Atom' object about the x axis
    # Arguments: (1) 'Atom' object
    #            (2) Angle by which to rotate, in radians
    # Modules:   numpy as np
    # Functions: NONE
    def rot_x(self, thetax=0):
        x = np.around(self.x, decimals=5)
        y = np.around((np.cos(thetax) * self.y) + (-1 * (np.sin(thetax) * self.z)), decimals=5)
        z = np.around((np.sin(thetax) * self.y) + (np.cos(thetax) * self.z), decimals=5)
        new_pos = Atom(x, y, z, self.mass, self.name, self.resn)
        return new_pos

    # Purpose:   To rotate an 'Atom' object about the y axis
    # Arguments: (1) 'Atom' object
    #            (2) Angle by which to rotate, in radians
    # Modules:   numpy as np
    # Functions: NONE
    def rot_y(self, thetay=0):
        x = np.around((np.cos(thetay) * self.x) + (np.sin(thetay) * self.z), decimals=5)
        y = np.around(self.y, decimals=5)
        z = np.around((np.sin(thetay) * self.x * -1) + (np.cos(thetay) * self.z), decimals=5)
        new_pos = Atom(x, y, z, self.mass, self.name, self.resn)
        return new_pos

    # Purpose:   To rotate an 'Atom' object about the z axis
    # Arguments: (1) 'Atom' object
    #            (2) Angle by which to rotate, in radians
    # Modules:   numpy as np
    # Functions: NONE
    def rot_z(self, thetaz=0):
        x = np.around((np.cos(thetaz) * self.x) - (np.sin(thetaz) * self.y), decimals=5)
        y = np.around((np.sin(thetaz) * self.x) + (np.cos(thetaz) * self.y), decimals=5)
        z = np.around(self.z, decimals=5)
        new_pos = Atom(x, y, z, self.mass, self.name, self.resn)
        return new_pos

    # Purpose:   To rotate an 'Atom' object about all three of the
    #            principal spatial axes by a specified amount
    # Arguments: (1) 'Atom' object
    #            (2) List object containing the values by which
    #                the x, y, and, z axes should be rotated about.
    #                must occur sequentially in the list; x, y, z
    #                order
    # Modules:   copy
    # Functions: 'rot_x' 'rot_y' 'rot_z'
    def rotation(self, theta=[0, 0, 0]):
        new_atom = copy.deepcopy(self)
        new_atom = new_atom.rot_x(thetax=theta[0])
        new_atom = new_atom.rot_y(thetay=theta[1])
        new_atom = new_atom.rot_z(thetaz=theta[2])
        return new_atom


'''----------
| Functions |
----------'''


# Purpose:   To read an xyz file into a dictionary format with all
#            of the atoms represent using 'Atom()' objects.
# Arguments: (1) File, with its absolute path included, as a string
# Modules:   string
# Functions: NONE
def read_xyz(file):
    with open(file, 'r') as f:
        contents = f.readlines()

    molecule = {}
    wspace = string.whitespace
    resn_name = contents[1].strip(wspace)
    contents = contents[2:]
    i = 1
    for line in contents:
        line = line.split()
        for j in range(len(line)):
            line[j] = line[j].strip(wspace)
        molecule[i] = Atom(float(line[1]), float(line[2]), float(line[3]),
                           atomic_masses[line[0][0]], line[0], resn_name)
        i += 1
    return molecule


# Purpose:   To read a pdb file into a dictionary format with all
#            of the atoms represent using 'Atom()' objects.
# Arguments: (1) File, with its absolute path included, as a string
# Modules:   string
# Functions: NONE
def read_pdb(file):
    with open(file, 'r') as f:
        contents = f.readlines()

    molecule = {}
    wspace = string.whitespace
    i = 1
    for line in contents:
        line = line.split()
        if line[0] == 'ATOM':
            for j in range(len(line)):
                line[j] = line[j].strip(wspace)
            molecule[i] = Atom(round(float(line[5]), 3),
                               round(float(line[6]), 3),
                               round(float(line[7]), 3),
                               atomic_masses[line[2][0]], line[2],
                               line[3])
            i += 1
        else:
            pass
    return molecule

# Purpose:   Calculates and returns the centroid of a molecule
# Arguments: (1) Dictionary representation of a molecule returned
#                by the 'read_xyz' function.
# Modules:   NONE
# Functions: NONE
def centroid(molecule):
    x, y, z = 0, 0, 0
    keys = list(molecule.keys())
    for atom in keys:
        x += molecule[atom].x
        y += molecule[atom].y
        z += molecule[atom].z
    div = len(keys)
    centroid = [x / div, y / div, z / div]
    return centroid


# Purpose:   To move the position of an entire molecule by a
#            specified about
# Arguments: (1) Dictionary representation of a molecule returned
#                by the 'read_xyz' function.
#            (2) A list containing the amount by which to shift the
#                x, y, and z coordinates.
# Modules:   NONE
# Functions: 'translate_atom'
def translate_molecule(molecule, t=[0, 0, 0]):
    new_molecule = copy.deepcopy(molecule)
    keys = list(new_molecule.keys())
    for atom in keys:
        new_molecule[atom] = new_molecule[atom].translate_atom(t)
    return new_molecule


# Purpose:   Returns the pair of H and O atoms with the greatest
#            distance between them. They will serve to define the
#            principle axis of the lipid. Returned object is a list
#            [index1, name1, index2, name2]
# Arguments: (1) A dictionary representation of the molecular
#                structure provided by the 'read_xyz' format
# Modules:   NONE
# Functions: NONE
def lipid_extrema(molecule):
    keys = list(molecule.keys())
    o_and_h = []
    for atom in keys:
        ele = molecule[atom].name[0]
        if ele == 'O' or ele == 'H':
            o_and_h.append(atom)
        else:
            pass
    i = 1
    distances = []
    for a1 in o_and_h:
        # The index prevents redundant comparisons
        for a2 in o_and_h[i:]:
            distances.append([a1, a2])
            i += 1
    for pair in distances:
        m1, m2 = pair[0], pair[1]
        dist = molecule[m1].distance(molecule[m2])
        pair.append(dist)
    distances.sort(key=lambda x: x[2])
    extrema = distances[-1][0:2]
    extrema = [extrema[0], molecule[extrema[0]].name, extrema[1], molecule[extrema[1]].name]
    if extrema[1][0] == 'O':
        extrema = [extrema[2], extrema[3], extrema[0], extrema[1]]
    else:
        pass
    return extrema


# Purpose:   To calculate the angle between two vectors.
# Arguments: (1) First vector, list or tuple object
#            (2) Second vector, list or tuple object
#            (3) Units which the angle value are returned in. The
#                default is 'rad' for radians and the other option
#                is 'deg' for degrees
# Modules:   numpy as np
# Functions: NONE
def vec_angle(v1, v2, units='rad'):
    angle = np.dot(np.array(v1), np.array(v2))
    angle = angle / (np.linalg.norm(np.array(v1)) * np.linalg.norm(np.array(v2)))
    angle = np.arccos(angle)
    if units == 'rad':
        pass
    elif units == 'deg':
        angle = (180 / np.pi)
    return angle


# Purpose:   To rotate a dictionary-representation of a molecule
#            about one of the principal axes, x, y, or z.
# Arguments: (1) A dictionary-representation of a molecule like that
#                returned by the 'read_xyz' function
#            (2) A list containing the values by which to rotate the
#                molecule about the principle axes. Values occur
#                sequentially, x, y, z.
# Modules:   copy
# Functions: 'rotation'
# NOTE:      The function does not accommodate multiple rotations.
#            You should only perform a rotation about one axis at a
#            time
def rot_molecule(molecule, rot=[0, 0, 0]):
    new_molecule = copy.deepcopy(molecule)
    keys = list(new_molecule.keys())
    for atom in keys:
        new_molecule[atom] = new_molecule[atom].rotation(theta=rot)
    return new_molecule


# Purpose:   To return the most extreme atom along a particular axis.
#            the extreme can be either high or low.
# Arguments: (1) A dictionary representation of a molecule such as
#                that provided by the 'read_xyz' function
#            (2) The axis along which to select the most extreme
#                atom; possible arguments are 'x' 'y' and 'z' with
#                'x' being the default
#            (3) The extreme can be either the highest or the lowest;
#                'lowest' for the lowest extreme and 'highest' for
#                the highest extreme
# Modules:   NONE
# Functions: NONE
def get_extreme(molecule, cor='x', extreme='lowest'):
    keys = molecule.keys()
    cors_list = []
    for atom in keys:
        cors_list.append([atom, molecule[atom].x, molecule[atom].y, molecule[atom].z])
    if cor == 'x':
        cors_list.sort(key=lambda x: x[1])
    elif cor == 'y':
        cors_list.sort(key=lambda x: x[2])
    elif cor == 'z':
        cors_list.sort(key=lambda x: x[3])

    if extreme == 'lowest':
        extreme_key = cors_list[0][0]
    elif extreme == 'highest':
        extreme_key = cors_list[-1][0]
    return extreme_key


# Purpose:   To write a dictionary-representation of a molecule,
#            such as that returned by the 'read_xyz' function as
#            a PDB file (God forgive me for this absolute
#            train wreck).
# Arguments: (1) A dictionary representation of a molecule such as
#                that returned by the 'read_xyz'
#            (2) The residue number of the molecule as a string
#                default is '1'; important for system formatting
#            (3) The starting atomic serial number;
# Modules:   NONE
# Functions: NONE
# NOTE:      This is a cluster-fuck. You really need to revisit this
#            function and simplify it, Samuel. Don't blame having to
#            work within someone else's jank-ass legacy format. Even
#            though a space-deliminated format would have been
#            vastly superior
def format_molecule_pdb(molecule, resn_seq='1', atom_serial='1'):
    keys = list(molecule.keys())
    keys.sort()
    record_name = 'ATOM'
    alternate_loc = ' '
    chain = ' '
    inser_code = ' '
    occ = '1.00'
    tmp_fac = '0.00'
    seg_iden = ''
    charge = ''
    all_entries = []
    for atom in keys:
        atom_name = molecule[atom].name
        resn_name = molecule[atom].resn
        x = molecule[atom].x
        y = molecule[atom].y
        z = molecule[atom].z
        elem = molecule[atom].name[0]
        entry = []
        # Record name
        entry.append(record_name + ((6 - len(record_name)) * ' '))
        # Atom serial number
        entry.append(((5 - len(atom_serial)) * ' ') + atom_serial)
        # Atom name
        entry.append('  ' + atom_name + ((3 - len(atom_name)) * ' '))
        # Alternate location
        entry.append(alternate_loc)
        # Residue name
        entry.append(resn_name + ' ')
        # Chain identifier
        entry.append(chain)
        # Residue sequence number
        entry.append(((4 - len(resn_seq)) * ' ') + resn_seq)
        # Code for the insertion of residues
        entry.append(inser_code)
        # x-coordinate
        x = str(round(x, 3))
        entry.append('   ' + ((8 - len(x)) * ' ') + x)
        # y-coordinate
        y = str(round(y, 3))
        entry.append(((8 - len(y)) * ' ') + y)
        # z-coordinate
        z = str(round(z, 3))
        entry.append(((8 - len(z)) * ' ') + z)
        # Occupancy
        entry.append(((6 - len(occ)) * ' ') + occ)
        # Temperature factor
        entry.append(((6 - len(tmp_fac)) * ' ') + tmp_fac)
        # Segment identifier
        entry.append('      ' + seg_iden + ((4 - len(seg_iden)) * ' '))
        # Elemental symbol
        entry.append(((2 - len(elem)) * ' ') + elem)
        # MAY NEED TO AMEND THIS GUY
        # Charge
        #
        entry.append(((2 - len(charge)) * ' ') + charge)
        entry.append('\n')
        entry = ''.join(entry)
        all_entries.append(entry)
        if int(atom_serial) >= 99999:
            pass
        else:
            atom_serial = str(int(atom_serial) + 1)
    return all_entries


# Purpose:
# Arguments:
# Modules:
# Functions:
def format_system_pdb(system, name):
    full_system = []
    system_keys = list(system.keys())
    resn = '1'
    serial = '1'
    for molecule in system_keys:
        molecule = system[molecule]
        molecule = format_molecule_pdb(molecule, resn_seq=resn, atom_serial=serial)
        ter_raw = molecule[-1].split()
        serial = str(int(ter_raw[1]) + 1)
        ter_line = ['TER   ']
        ter_line.append(((5 - len(serial))) * ' ' + serial)
        ter_line.append('      ' + ter_raw[3])
        ter_line.append('  ')
        ter_line.append(((4 - len(ter_raw[4])) * ' ') + ter_raw[4])
        ter_line.append('\n')
        ter_line = ''.join(ter_line)
        molecule.append(ter_line)
        full_system.append(molecule)
        serial = str(int(serial) + 1)
        resn = str(int(resn) + 1)
    full_system.append(['END\n'])
    #return full_system


    # Purpose:
    # Arguments:
    # Modules:
    # Functions:
    #def write_formatted_system(formatted_system, name)
    with open(name, 'w') as f:
        for molecule in full_system:
            for atom in molecule:
                f.writelines(atom)



'''-----------
| Main Block |
-----------'''

def micelle_scaffold()

# Import file contents
poc_file = '/home/sam/Desktop/model_proj/programs/MMAEVe/lipids/POC.pdb'
lipid = read_pdb(poc_file)

# Move first extrema to the origin
extreme_points = lipid_extrema(lipid)
e1 = copy.deepcopy(lipid[extreme_points[0]])
origin = Atom()
trans = origin - e1
lipid = translate_molecule(lipid, t=trans)

# Update extrema values
e1 = copy.deepcopy(lipid[extreme_points[0]])
e2 = copy.deepcopy(lipid[extreme_points[2]])

# Move second extrema to the x-axis
# The angles are multiplied by a correction factor that forces them
# to be rotated in the proper direction. Correction is achieved by
# checking a sign.
z_angle = vec_angle([e2.x, e2.y, 0], [1, 0, 0]) * np.sign(e2.y) * -1
y_angle = vec_angle([e2.x, 0, e2.z], [1, 0, 0]) * np.sign(e2.z) * -1
rot_angles = [0, 0, z_angle]
lipid = rot_molecule(lipid, rot=rot_angles)
rot_angles = [0, y_angle, 0]
lipid = rot_molecule(lipid, rot=rot_angles)

# Update extrema values
e1 = copy.deepcopy(lipid[extreme_points[0]])
e2 = copy.deepcopy(lipid[extreme_points[2]])

# Shift the aligned molecule so that it is to the right of the
# y-axis
low_x = get_extreme(lipid, cor='x', extreme='lowest')
low_x = lipid[low_x].x
trans = 0 - low_x
lipid = translate_molecule(lipid, t=[trans + 1, 0, 0])

# Create the initial semi-circle of the micelle
semi_circle = {}
radius = get_extreme(lipid, cor='x', extreme='highest')
radius = lipid[radius].x
semi_circumference = radius * np.pi
sections = int(semi_circumference / 8)
sec_length = semi_circumference / sections
sec_rot = sec_length * (np.pi / semi_circumference)

i = 1
semi_circle[i] = lipid
# Sections + 2 to add an extra lipid to the initial semi circle.
# This ensures that there are no significant gaps in the initial
# structure
for i in range(2, sections + 2):
    semi_circle[i] = rot_molecule(semi_circle[i - 1], rot=[0, 0, sec_rot])

# Create the full circle using the semi-circle
# TEST: Currently testing to see if the high_y will work for
#       defining the circumference of each slice
#
full_circle = {}
semi_circle_keys = list(semi_circle.keys())
semi_circle_keys.sort()
i = 2
for slice in semi_circle_keys:
    if slice == 1:
        sec_div = 24
    elif slice == 2:
        sec_div = 48
    elif slice == 3:
        sec_div = 16
    else:
        sec_div = 10
    high_y = get_extreme(semi_circle[slice], cor='y', extreme='highest')
    radius = semi_circle[slice][high_y].y
    circumference = radius * 2 * np.pi
    sections = int(circumference / sec_div)
    sec_length = circumference / sections
    sec_rot = sec_length * ((2 * np.pi) / circumference)

    full_circle[i] = semi_circle[slice]
    i += 1
    k = i
    for j in range(k, k + sections -1):
        full_circle[j] = rot_molecule(full_circle[i - 1], rot=[sec_rot, 0, 0])
        i += 1

# Test
format_system_pdb(full_circle, name='/home/sam/Desktop/model_proj/programs/MMAEVe/lipids/test.pdb')
