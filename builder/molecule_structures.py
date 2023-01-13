'''------------------------------------------------------*
| Title:   Molecule Structures                           |
|                                                        |
| Author:  Gubbin Eel (Satanic Overlord of the Swamp)    |
|                                                        |
| Purpose: Defines classes needed to represent molecules |
|          along with the methods needed to effectively  |
|          manipulate them.                              |
*------------------------------------------------------'''

import geom_shapes as gs
import copy
import numpy as np

# Atomic masses of common elements.
atomic_masses = {'H' : 1.00794, 'C' : 12.0107, 'O' : 15.9994, 
                 'N' : 14.0067, 'F' : 18.9984, 'P' : 30.9738, 
                 'S' : 32.0650, 'Cl': 35.4530, 'Br': 79.9040, 
                 'I' : 126.904, 'He': 4.00260, 'Ne': 20.1797, 
                 'Ar': 39.9480, 'Li': 6.94100, 'Be': 9.01218, 
                 'B' : 10.8110, 'Na': 22.9898, 'Mg': 24.3050, 
                 'Al': 26.9815, 'Si': 28.0855, 'K' : 39.0983, 
                 'Ca': 40.0780, 'Sc': 44.9559, 'Ti': 47.8670,
                 'V' : 50.9415, 'Cr': 51.9961, 'Mn': 54.9380, 
                 'Fe': 55.8450, 'Co': 58.9332, 'Ni': 58.6934, 
                 'Cu': 63.5460, 'Zn': 65.4090, 'Ga': 69.7230, 
                 'Ge': 72.6400, 'As': 74.9216, 'Se': 78.9600, 
                 'Kr': 83.7980, 'X' : 0.00000, ''  : 72.0000}
# Atomic radii of common elements.
atomic_radii = {  'H' : 0.37, 'C' : 0.77, 'O' : 0.73, 'N' : 0.75, 
                  'F' : 0.71, 'P' : 1.10, 'S' : 1.03, 'Cl': 0.99, 
                  'Br': 1.14, 'I' : 1.33, 'He': 0.30, 'Ne': 0.84, 
                  'Ar': 1.00, 'Li': 1.02, 'Be': 0.27, 'B' : 0.88, 
                  'Na': 1.02, 'Mg': 0.72, 'Al': 1.30, 'Si': 1.18, 
                  'K' : 1.38, 'Ca': 1.00, 'Sc': 0.75, 'Ti': 0.86, 
                  'V' : 0.79, 'Cr': 0.73, 'Mn': 0.67, 'Fe': 0.61, 
                  'Co': 0.64, 'Ni': 0.55, 'Cu': 0.46, 'Zn': 0.60, 
                  'Ga': 1.22, 'Ge': 1.22, 'As': 1.22, 'Se': 1.17, 
                  'Kr': 1.03, 'X' : 0.00, ''  : 2.35}

class Atom(object):
    '''
    Description: 
        Represent a single atom that is intended to be a member of a 
        more complex structure or molecule. Contains important 
        information needed to describe atoms. As well as the basics 
        for writing the major molecular file-formats.

    Attributes:  
        serial) Integer. Atom serial number. Common field in most 
        popular biomolecule formats.
        x) Float. x-coordinate.
        y) Float. y-coordinate.
        z) Float. z-coordinate.
        mass) Float. Atomic mass.
        elem) String. Element name.
        name) String. Atom name.
        resn) String. Residue name. 
        resi) Integer. Residue number.
    '''

    def __init__(self,        serial = 0,  x = 0.,      y = 0., 
                 z = 0.,      mass = 0.,   resi = 0,    resn = "NA", 
                 chain = " ", elem = "NA", name = "NA", charge = 0.):
        ''' 
        Purpose:   Initialize an Atom instance.
        Arguments: self) Atom instance.
                   serial) Integer. Atom serial number. Common field 
                   in most populat biomolecule formats.
                   x) Float. x-coordinate.
                   y) Float. y-coordinate.
                   z) Float. z-coordinate.
                   mass) Float. Atomic mass.
                   resi) Integer. Residue number.
                   resn) String. Residue name.
                   chain) String. Chain biomolecule belongs to.
                   elem) String. Elemental symbol of the atom.
                   name) String. Atom name of residue.
                   charge) Float. Charge on the atom.
        Returns:   Atom instance.
        '''
        self.serial = serial
        self.x      = x   
        self.y      = y   
        self.z      = z   
        self.mass   = mass
        self.resi   = resi
        self.resn   = resn
        self.chain  = chain
        self.elem   = elem
        self.name   = name
        self.charge = charge

    def translate(self, t = [0., 0., 0,]):
        ''' 
        Purpose:   Translate Atom in the x, y, and z-directions  
                   using a provided vector.
        Arguments: 1) Atom instance.
                   2) List or other iterable. Should contain the 
                      x, y, and z values, respectively, that the
                      the molecule will be translated by.
        Requires:  Nothing
        Modifies:  x, y, z
        Returns:   Nothing
        '''
        self.x += t[0]
        self.y += t[1]
        self.z += t[2]

    def distance(self, other):
        ''' 
        Purpose:   Calculates the distance between two atoms.
        Arguments: self) Atom instance.
                   other) Second instance of the Atom class.
        Requires:  x, y, z
        Modifies:  Nothing
        Returns:   Float. The distance between two atoms.
        '''
        dist = (((self.x - other.x) ** 2)  +          \
                ((self.y - other.y) ** 2)  +          \
                ((self.z - other.z) ** 2)) ** (1 / 2)
        return dist

    def distance_vec_atom(self, other):
        ''' 
        Purpose:   Generates a vector in the form of a list that 
                   points toward another atom.
        Arguments: 1) Atom instance.
                   2) Second instance of the Atom class.
        Requires:  x, y, z
        Modifies:  Nothing
        Returns:   List of Floats. Contains the x, y, and z 
                   coordinates, respectively of a distance vector 
                   pointing from the position of an Atom instance to 
                   the position of another Atom instance.
        '''
        dist_vec = [other.x - self.x,
                    other.y - self.y,
                    other.z - self.z]
        return dist_vec

    def distance_vec_point(self, point):
        ''' 
        Purpose:   Generates a vector in the form of a list that 
                   points toward a point.
        Arguments: self) Atom instance.
                   point) 1 x 3 List or other iterable of Floats. 
                   Contains x, y, and z coordinates, respectively, of 
                   the point to be translated to. 
        Requires:  x, y, z
        Modifies:  Nothing
        Returns:   List of Floats. Contains the x, y, and z 
                   coordinates, respectively, of a distance vector 
                   pointing from the current position of the Atom 
                   instance to the specified point.
        '''
        dist_vec = [point[0] - self.x,
                    point[1] - self.y,
                    point[2] - self.z]
        return dist_vec

    def rot_x(self, thetax = 0.):
        ''' 
        Purpose:   Rotate an Atom object about the x-axis.
        Arguments: self) Atom instance.
                   thetax) Float. Angle by which to rotate the atom, 
                   in radians.
        Requires:  y, z
        Modifies:  y, z
        Returns:   None
        '''
        y = self.y
        z = self.z
        self.y = (np.cos(thetax) * y) + (-1.0 * (np.sin(thetax) * z))
        self.z = (np.sin(thetax) * y) + (np.cos(thetax) * z)

    def rot_y(self, thetay = 0.):
        ''' 
        Purpose:   Rotate an Atom object about the y-axis.
        Arguments: self) Atom instance.
                   thetay) Float. Angle by which to rotate the atom, 
                   in radians.
        Requires:  x, z
        Modifies:  x, z
        Returns:   None
        '''
        x = self.x
        z = self.z
        self.x = (np.cos(thetay) * x) + (np.sin(thetay) * z)
        self.z = (np.sin(thetay) * x * -1.0) + (np.cos(thetay) * z)

    def rot_z(self, thetaz = 0.):
        ''' 
        Purpose:   Rotate an Atom object about the z-axis.
        Arguments: self) Atom instance.
                   thetay) Float. Angle by which to rotate the atom, 
                   in radians.
        Requires:  x, y
        Modifies:  x, y
        Returns:   None
        '''
        x = self.x
        y = self.y
        self.x = (np.cos(thetaz) * x) - (np.sin(thetaz) * y) 
        self.y = (np.sin(thetaz) * x) + (np.cos(thetaz) * y)

    def to_pos_vec(self):
        ''' 
        Purpose:   Retrieve the coordinates of an Atom object as a 
                   list.
        Arguments: self) Atom instance.
        Requires:  x, y, z
        Modifies:  Nothing
        Returns:   List of Floats. Contains the x, y, and z 
                   coordinates, respectively, of an atom.
        '''
        return [self.x, self.y, self.z]

    def to_pdb_file(self, atom_serial = '0', resn_seq = '0', 
                    chain = ''):
        ''' 
        Purpose:   Formats an Atom instance as a PDB file line.
        Arguments: self) Atom instance.
                   atom_serial) Integer. The atom serial number. If it 
                   is '0' then the value already assigned will be 
                   used.
                   resn_seq) Integer. The residue sequence number. If 
                   it is 0 the value already assigned will be used.
                   chain) String. Chain ID. Default is nothing.
        Requires:  serial, name, resn, resi, chain, x, y, z, elem, 
                   charge
        Modifies:  Nothing
        Returns:   String. Formatted with all atom info specified in 
                   the proper .pdb format for an "ATOM" line.
        '''
        # All the required entries for an ATOM field entry in a pdb 
        # file.
        record_name = "ATOM"
        if atom_serial == '0':
            atom_serial = str(self.serial)
        else:
            pass
        if int(atom_serial) > 99999:
            atom_serial = "99999"
        atom_name = self.name
        alternate_loc = ' '
        resn_name = self.resn 
        chain = self.chain
        if resn_seq == '0':
            resn_seq = str(self.resi)
        else:
            pass
        inser_code = ' '
        x = str(format(round(self.x, 3), ".3f"))
        y = str(format(round(self.y, 3), ".3f"))
        z = str(format(round(self.z, 3), ".3f"))
        occ = '1.00'
        tmp_fac = '0.00'
        seg_iden = ' '
        elem = self.elem
        charge = str(self.charge)
        line = []
        # Record name
        line.append(record_name + ((6 - len(record_name)) * ' '))
        # Atom serial number
        line.append(((5 - len(atom_serial)) * ' ') + atom_serial)
        # Atom name
        line.append(' ' + atom_name + ((4 - len(atom_name)) * ' ')) 
        # Alternate location
        line.append(alternate_loc)
        # Residue name
        line.append(resn_name + ' ')
        # Chain identifier
        line.append(chain)
        # Residue sequence number
        line.append(((4 - len(resn_seq)) * ' ') + resn_seq)
        # Code for the insertion of residues
        line.append(inser_code)
        # x-coordinate
        line.append('   ' + ((8 - len(x)) * ' ') + x)
        # y-coordinate
        line.append(((8 - len(y)) * ' ') + y)
        # z-coordinate
        line.append(((8 - len(z)) * ' ') + z)
        # Occupancy
        line.append(((6 - len(occ)) * ' ') + occ)
        # Temperature factor
        line.append(((6 - len(tmp_fac)) * ' ') + tmp_fac)
        # Segment identifier
        line.append('   ' + seg_iden + ((4 - len(seg_iden)) * ' '))
        # Elemental symbol
        line.append(((2 - len(elem)) * '    ') + elem)
        # Charge
        line.append(((2 - len(charge)) * ' ') + charge)
        line.append('\n')
        line = ''.join(line)
        return line

class Molecule(object):
    '''
    Description: 
        Contains the constituent atoms of single molecule along with 
        relevant information for operations, primarily writing various 
        file formats and performing spatial manipulations on the 
        molecular coordinates.

    Attributes:  
        serial) N x 1 NP Array of Integers. Serial number of the atom.
        xyz) N x 3 NP Array of Floats. x, y, and z coordinates of the 
        atoms, respectively.
        mass) N x 1 NP Array of Floats. Atomic masses of the atoms.
        resi) N x 1 NP Array of Integers. The residue numbers of the 
        atoms.
        resn) N x 1 NP Array of Strings. The residue name of the 
        atoms.
        chain) N x 1 NP Array of Strings. Contains the chain IDs of 
        the atoms.
        elem) N x 1 NP Array of Strings. Contains the elemental 
        symbols of the atoms.
        name) N x 1 NP Array of Strings. Contains the names of the 
        atoms. 
        charge) N x 1 NP Array of Floats. Contains the charge on 
        each atom.
        radius) N x 1 NP Array of Floats. Contains the atomic radii of 
        the atoms.
        lookup) Dictionary. Uses unique keys for the quick lookup of 
        an atom of interest. The keys are built using the atom serial 
        numbers, residue numbers, and residue names seperated by '-'.
        extrema) Dictionary. Meant to contain the indicies of the 
        information for the polar extrema of the molecule. Intended 
        use is by the Lipid child class. It will store which atom is 
        the "Tail" and which is the "Head" of the Lipid/Molecule. 
    '''

    def __init__(self,                   serial = np.array([]), 
                 xyz     = np.array([]), mass   = np.array([]),
                 resi    = np.array([]), resn   = np.array([]),   
                 chain   = np.array([]), elem   = np.array([]),   
                 name    = np.array([]), charge = np.array([]), 
                 radius  = np.array([]), lookup = dict(),
                 extrema = dict()):
        ''' 
        Purpose:   To initialize a Molecule instance.
        Arguments: self) Molecule instance.
                   serial) 1 x N NP array of integers. The serial 
                   numbers of the individual atoms.
                   xyz) 3 x N NP array of floats. The x, y, and z 
                   coordinates, respectively, of the atoms.
                   mass) 1 x N NP array of floats. The masses of the 
                   atoms.
                   resi) 1 x N NP array of integers. The residue 
                   numbers of the atoms, (usually all the same within) 
                   a given smallish molecule.
                   resn) 1 x N NP array of strings. The residue name.
                   chain) 1 x N NP array of strings. The chain id.
                   elem) 1 x N NP array of strings. The element name.
                   name) 1 x N NP array of strings. The atom names.
                   charge) 1 x N NP array of floats. Atomic charge.
                   radius) 1 x N NP array of floats. Atomic radii.
                   lookup) Dictionary. Dictionary with the key as a 
                   combination of the residue name, number, and atom 
                   number.                   
                   extrema) Dictionary. Dictionary with two keys 
                   "Head" and "Tail". These provide the index of the 
                   atoms which serve correspond to the orientation of 
                   a Lipid's head and tail in a biological membrane.
        Requires:  Nothing
        Modifies:  All
        Returns:   Molecule instance.
        '''
        self.serial  = serial
        self.xyz     = xyz
        self.mass    = mass
        self.resi    = resi
        self.resn    = resn
        self.chain   = chain
        self.elem    = elem
        self.name    = name
        self.charge  = charge
        self.radius  = radius
        self.lookup  = lookup
        self.extrema = extrema

    def translate(self, t = np.array([0., 0., 0.])):
        ''' 
        Purpose:   Translates a molecule in the x, y, and z directions
                   using a provided list.
        Arguments: self) Molecule instance.
                   t) 1 x 3 NP Array or other iterable. Contains 
                   coordinates of the x, y, and z positions, 
                   respectively. Describes a vector used to translate 
                   the molecule.
        Requires:  Nothing
        Modifies:  xyz
        Returns:   Nothing. Updates Molecule instance.
        '''
        self.xyz += t

    def rotate(self, axis = 'x', theta = 0.):
        ''' 
        Purpose:   To rotate a molecule by a specified amount around 
                   either the x, y, or z-axis.
        Arguments: self) Molecule instance.
                   axis) String. The axis to rotate about: 'x', 'y', 
                   or 'z'. 
                   theta) Float. How far to rotate, in radians; 
                   default of 0.
        Requires:  xyz
        Modifies:  xyz
        Returns:   None. Updates Molecule instance.
        '''
        if axis == 'x':
            self.xyz[:, 1],                                     \
            self.xyz[:, 2] =                                    \
                             (np.cos(theta) * self.xyz[:, 1]) - \
                             (np.sin(theta) * self.xyz[:, 2])   \
                             ,                                  \
                             (np.sin(theta) * self.xyz[:, 1]) + \
                             (np.cos(theta) * self.xyz[:, 2])
        if axis == 'y':
            self.xyz[:, 0],                                     \
            self.xyz[:, 2] =                                    \
                             (np.cos(theta) * self.xyz[:, 0]) + \
                             (np.sin(theta) * self.xyz[:, 2])   \
                             ,                                  \
                             (np.cos(theta) * self.xyz[:, 2]) - \
                             (np.sin(theta) * self.xyz[:, 0])   
        if axis == 'z':
            self.xyz[:, 0],                                     \
            self.xyz[:, 1] =                                    \
                             (np.cos(theta) * self.xyz[:, 0]) - \
                             (np.sin(theta) * self.xyz[:, 1])   \
                             ,                                  \
                             (np.sin(theta) * self.xyz[:, 0]) + \
                             (np.cos(theta) * self.xyz[:, 1])   

    def to_pdb_lines(self, atom_start = 1, res_num = 1):
        '''
        Purpose:   Converts a molecule to a list of strings formatted
                   for a pdb file.
        Arguments: self) Molecule instance.
                   atom_start) The starting number for numbering 
                   atoms.
                   res_num) The residue number.
        Requires:  serial, name, resn, chain, xyz, elem, charge
        Modifies:  Nothing
        Returns:   List of Strings. Contains strings for each of the 
                   atoms of the molecule formatted as pdb lines.
        '''
        record_name   = np.array(["ATOM" for ii in self.serial])
        serial = np.arange(atom_start, 
                           atom_start + len(self.serial))
        serial[serial > 99999] = 99999
        serial        = serial.astype(str)
        atom_name     = self.name
        alternate_loc = np.array([' ' for ii in self.serial])
        residue_name  = self.resn
        chain         = self.chain
        residue_num   = np.ones(len(self.resn)) * res_num
        residue_num[res_num > 9999] = 9999
        residue_num   = residue_num.astype(int).astype(str)
        insertion     = alternate_loc
        x             = np.array(["%.3f" % round(xx, 3) \
                                  for xx in self.xyz[:, 0]])
        y             = np.array(["%.3f" % round(yy, 3) \
                                  for yy in self.xyz[:, 1]])
        z             = np.array(["%.3f" % round(zz, 3) \
                                 for zz in self.xyz[:, 2]])
        occupancy     = np.array(['1.00' for ii in self.serial])
        temp_fac      = np.array(['0.00' for ii in self.serial])
        seg_iden      = alternate_loc
        element       = self.elem
        charge        = self.charge.astype(str)
        dummy_field   = np.array(['' for ii in self.serial])

        all_fields = np.array([record_name,  dummy_field, 
                               serial,       dummy_field,  
                               atom_name,    alternate_loc,
                               residue_name,  
                               chain,        residue_num,  
                               insertion,    dummy_field,
                               x,            y,           
                               z,            occupancy,    
                               temp_fac,     dummy_field,
                               seg_iden,     element,     
                               charge])
        field_size = [4, 2, 5, 1, 4, 1, 4, 1, 4, 1, 3, 8, 8,
                      8, 6, 6, 6, 4, 2, 2]
        side = ['left',  'left',  'right', 'left',  'left',  'left', 
                'left',  'left',  'right', 'left',  'left',
                'right', 'right', 'right', 'right', 'right', 'left',
                'left',  'left',  'right']
        
        table_space = alternate_loc
        for ii in range(len(all_fields)):
            field_lengths = \
                    np.array([len(entry) for entry in all_fields[ii]])
            space_length = field_size[ii] - field_lengths
            space_to_add = np.char.multiply(' ', space_length)
            if side[ii] == 'left':
                all_fields[ii] = np.char.add(all_fields[ii], 
                                             space_to_add)
            else:
                all_fields[ii] = np.char.add(space_to_add, 
                                             all_fields[ii])
        pdb_lines = np.apply_along_axis(''.join, 0, all_fields)
        pdb_lines = np.char.add(pdb_lines, '\n')
        pdb_lines = pdb_lines.tolist()

        return(pdb_lines)

    def add_mass(self):
        '''
        Purpose:   Add atomic masses to Molecule instance.
        Arguments: self) Molecule instance.
        Requires:  elem
        Modifies:  mass
        Returns:   Nothing
        '''
        self.mass = np.array([atomic_masses[elem]    \
                             for elem in self.elem])

    def add_radius(self):
        '''
        Purpose:   Add atomic radii to Molecule instance.
        Arguments: self) Molecule instance.
        Requires:  elem
        Modifies:  radius
        Returns:   None. Updates Molecule instance.
        '''
        #try:
        #    self.radius = np.array([atomic_radii[elem]      \
        #                           for elem in self.elem ])
        #except KeyError:
        #    self.radius = np.array([atomic_radii[name]      \
        #                           for name in self.name])
        self.radius = np.array([atomic_radii[elem]      \
                               for elem in self.elem ])
        
    def generate_lookup(self):
        '''
        Purpose:   Generate/update the dictionary for looking up an 
                   atom's index.
        Arguments: self) Molecule instance.
        Requires:  serial, resi, resn
        Modifies:  lookup
        Returns:   Nothing
        '''
        for ii in range(len(self.serial)):
            key = str(int(self.serial[ii])) + '-' + \
                  str(int(self.resi[ii]))   + '-' + \
                      self.resn[ii]
            self.lookup[key] = ii

class Lipid(Molecule):
    '''
    Description:
        Inherits Molecule class. Largely similar to the Molecule class.
        New methods are related to orienting Lipid instances within 
        a given structure.

    Attributes:
        See the Molecule class description.
    '''

    def lipid_extrema(self): 
        ''' 
        Purpose:   Determines which atoms/beads should be used as the 
                   head and tail ends of the Lipid. These atoms are 
                   used for aligning the Lipids when building the 
                   various structures. Information is stored as a 
                   Dictionary in the extrema attribute.
        Arguments: self Lipid instance.
        Requires:  name, xyz
        Modifies:  extrema
        Returns:   Nothing
        '''
        # List to hold Atom class instances for negative atoms/beads.
        head_atoms = [] 
        # List to hold Atom class instances for positive atoms/beads.
        tail_atoms = [] 
        # Retrieve class instances from the dictionary.
        ii = 0
        for name in self.name:
            if any(x in name for x in ['G', 'P', 'O', 'N']):
                head_atoms.append(ii)
            elif any(x in name for x in ['C', 'D', 'H']):
                tail_atoms.append(ii)
            else:
                pass
            ii += 1

        # List of lists to contain the atomic serial numbers and 
        # their distance.
        distance_comp = []
        head_len, head_len = len(head_atoms), len(tail_atoms)
        for head in head_atoms:
            for tail in tail_atoms:
                # Distance between h and o atoms.
                dist    = np.linalg.norm(self.xyz[head, :] - \
                                         self.xyz[tail, :])
                # Molecule dictionary key.
                head_index = head
                # Molecule dictionary key.
                tail_index = tail
                distance_comp.append([head_index, tail_index, dist])

        # Sort the list from largest to shortest distance.
        distance_comp = sorted(distance_comp, key = lambda x: x[2])
        distance_comp = distance_comp[-1]
        distance_comp = {"Head" : distance_comp[0], 
                         "Tail" : distance_comp[1], 
                         "dist" : distance_comp[2]}
        self.extrema = distance_comp

    def trans_axis_to_point(self, point, axis_member):
        ''' 
        Purpose:   Translate either the Lipid head of tail to a point.
        Arguments: self) Lipid instance.
                   point) 1 x 3 NP array of floats. Contains the x, y, 
                   and z coordinates, respectively, of the point that 
                   the Head or Tail should be translated to.
                   axis_member) String. The name of the axis member to 
                   translate to the piont, "Head" or "Tail".
        Requires:  extrema, xyz
        Modifies:  xyz  
        Returns:   Nothing
        '''
        # Get axis extreme, either Head or Tail.
        extreme = self.extrema[axis_member]
        # Distance vector required to translate oxygen to a point.
        distance_vec = point - self.xyz[extreme, :]
        # Translate the lipid to the position on the sphere.
        self.translate(t = distance_vec)

    def align_axis_to_vec(self, point, axis_member):
        ''' 
        Purpose:   To rotate one extrema member, either "Head" or 
                   "Tail", while holding the other fixed such that the 
                   member being rotated lands on the position vector 
                   specified by the point.
        Arguments: self) Lipid instance.
                   point) 1 x 3 NP Array. Position vector to align to.
                   axis_member) String. The name of the axis member to 
                   rotate, "Head" or "Tail".
        Requires:  xyz
        Modifies:  xyz
        Returns:   Nothing
        '''
        # Perform a rotation about each axis to align with the vector.
        # The rotation is performed by projecting along the rotation
        # axis and getting the angle between the two projections with
        # respect to the axis normal.

        # Select which axis extreme will be mobile. 
        extreme = self.extrema[axis_member]

        # z-axis rotation.
        sphere_vec = [point[0], point[1], 0]
        lipid_vec  = [self.xyz[extreme, 0], 
                      self.xyz[extreme, 1], 
                      0]
        angle = gs.directional_angle(lipid_vec, sphere_vec, [0, 0, 1])
        self.rotate(axis = 'z', theta = angle)

        # y-axis rotation.
        sphere_vec = [point[0], 0, point[2]]
        lipid_vec  = [self.xyz[extreme, 0], 
                      0, 
                      self.xyz[extreme, 2]]
        angle = gs.directional_angle(lipid_vec, sphere_vec, [0, 1, 0])
        self.rotate(axis = 'y', theta = angle)

        # x-axis rotation.
        sphere_vec = [0, point[1], point[2]]
        lipid_vec  = [0, 
                      self.xyz[extreme, 1], 
                      self.xyz[extreme, 2]]
        angle = gs.directional_angle(lipid_vec, sphere_vec, [1, 0, 0])
        self.rotate(axis = 'x', theta = angle)
 
        # Hello Darkness my old friend. You might ask yourself: 
        # "Weren't the rotations about the z and y axes already
        # performed?" You would be correct but, I encountered the
        # ancient nemesis of many a poor programmer, numerical error.
        # Turns out that using arccos has an associated rounding 
        # error, expecially when an angle is close to 0. So after
        # performing a rotation the result is that the other angles
        # get a little fucked. I spent a few hours looking into
        # alternatives then finally said: "FUCK IT" and tried 
        # applying the earlier two rotations again. Wouldn't you
        # know that it actually fucking worked. However, this isn't
        # an elegant solution by any stretch. I'm leaving this big
        # ass brick of a comment to draw attention to my (Gubbin) 
        # failure so that I remember to revisit it at a later 
        # date. Or maybe I get lucky and Dumpster Mokey takes care
        # of it. Though I hope for the former as the "Dumpster 
        # Monkey's Paw" phenomenon can be difficult to deal with; 
        # Dumpster Monkey will usually put out one fire only to 
        # light a dozen more in its place.

        # z-axis rotation.#
        sphere_vec = [point[0], point[1], 0]#
        lipid_vec  = [self.xyz[extreme, 0], #
                      self.xyz[extreme, 1], #
                      0]#
        angle = gs.directional_angle(lipid_vec, sphere_vec, [0, 0, 1])#
        self.rotate(axis = 'z', theta = angle)#
        # y-axis rotation.#
        sphere_vec = [point[0], 0, point[2]]#
        lipid_vec  = [self.xyz[extreme, 0], #
                      0, #
                      self.xyz[extreme, 2]]#
        angle = gs.directional_angle(lipid_vec, sphere_vec, [0, 1, 0])#
        self.rotate(axis = 'y', theta = angle)#
 
class Residue(Molecule):
    '''
    Description:
        Extension of the molecule class. Used to represent Protein 
        Residues in a syntactically distinct manner from Molecules 
        which is meant to represent plain old molecules.
    
    Attributes:
        See the description in the Molecule class.
    '''
    pass

class Protein(object):
    '''
    Description:
        Used to represent a proteins structure. Defines methods for 
        spatial manipulation and conversion to a .pdb file.

    Attributes: 
        residues) N x 1 NP Array of Residue instances. Contains the 
        residues of the protein as Residue instances. Otherwise, an 
        empty NP Array is assigned.
        extrema) Dictionary. Contains information about the extrema of 
        the protein structure. They are not, in-fact, extrema. They 
        are user-specified atoms. The naming convention is kept 
        because the user-specified "Head" and "Tail" are aligned with 
        the head and tail of the Lipid instances for each structure. 
        The key is either "Head" or "Tail" and the stored value is a 
        list or other iterable with two values. The first being the 
        residue index and the second being the atom index within the 
        residue.
    '''
    def __init__(self, residues = np.array([]), 
                 extrema = dict()):
        '''
        Purpose:   Initialize an instance of the Protein class.
        Arguments: self) Protein instance.
                   residues) N x 1 NP Array of Residue instances. 
                   Contains the residues of the protein as Residue 
                   instances. Otherwise, an empty NP Array is 
                   assigned.
                   extrema) Dictionary. Contains information about the  
                   extrema of the protein structure. They are not, 
                   in-fact, extrema. They are user-specified atoms. 
                   The naming convention is kept because the  
                   user-specified "Head" and "Tail" are aligned with 
                   the head and tail of the Lipid instances for each 
                   structure. The key is either "Head" or "Tail" and 
                   the stored value is a list or other iterable with 
                   two values. The first being the residue index and 
                   the second being the atom index within the residue.
        Returns:   Protein instance.
        '''
        self.residues = residues
        self.extrema  = extrema

    def translate(self, t = np.array([0., 0., 0.])):
        ''' 
        Purpose:   Translates a Protein instance in the x, y, and z 
                   directions using a provided vector.
        Arguments: self) Protein instance.
                   t) 1 x 3 NP Array. Contains three floats, x, y, 
                      and z, respectively, defining a direction
                      vector for translating the Protein instance.
        Requires:  residues
        Modifies:  residues
        Returns:   Nothing
        '''
        for res in self.residues:
            res.translate(t)

    def rotate(self, axis = 'x', theta = 0):
        ''' 
        Purpose:   To rotate a Protein by a specified amount around 
                   either the x, y, or z-axis.
        Arguments: self) Protein instance.
                   axis) String. The axis to rotate about: 'x', 'y', 
                   or 'z'. 
                   theta) Float. How far to rotate, in radians; 
                   default of 0.
        Requires:  residues
        Modifies:  residues
        Returns:   Nothing
        '''
        for res in self.residues:
            res.rotate(axis = axis, theta = theta)

    def protein_extrema(self, head_id, tail_id): 
        '''
        Purpose:   To extract the location information of the
                   specified atoms for quickly navigating
                   through the protein and assigning them to
                   extrema. The head and tail names are derived
                   from the fact that the selected atoms will
                   be aligned with the Lipid head and tail 
                   orientations, respectively.
        Arguments: self) Protein instance.
                   head_id) String. Head atom's serial number, residue
                   number, and residue name, seperated by
                   dashes. For example, '80054-911-DID', where 
                   '80054' is the serial number, '911' is the 
                   residue number, and 'DID' is the residue name.
                   tail_id) String. Tail atom's serial number, residue
                   number, and residue name, seperated by
                   dashed. For example, '420-69-VAG' where '420'
                   is the atom serial number, '69' is the 
                   residue number, and 'VAG' is the residue name.
        Requires:  residues
        Modifies:  extrema
        Returns:   Nothing
        '''
        res_number   = 0
        for res in self.residues:
            if head_id in res.lookup:
                self.extrema["Head"] = [res_number, 
                                        res.lookup[head_id]]
            elif tail_id in res.lookup:
                self.extrema["Tail"] = [res_number, 
                                        res.lookup[tail_id]]
            else:
                pass
            res_number += 1

    def trans_axis_to_point(self, point, axis_member):
        ''' 
        Purpose:   Translate either the tail or head atom of a Lipid 
                   to a point.
        Arguments: self) Protein instance.
                   point) 1 x 3 NP Array. Contains the x, y, and z 
                   coordinates, respectively, of the point that Tail 
                   atom should be translated to.
                   axis_member) String. The name of the axis member to 
                   rotate, "Head" or "Tail".
        Requires:  extrema, residues
        Modifies:  residues
        Return:    None. Updates Protein instance.
        '''
        # Get axis extreme, either Head or Tail.
        res  = self.extrema[axis_member][0]
        atom = self.extrema[axis_member][1]
        # Distance vector required to translate a given atom to a 
        # point.
        distance_vec = point - self.residues[res].xyz[atom, :]
        # Translate the lipid along a provided vector.
        self.translate(t = distance_vec)

    def align_axis_to_vec(self, point, axis_member):
        ''' 
        Purpose:   Align a Protein with a given position vector with 
                   either the tail or head atom remaining in its 
                   current position. 
        Arguments: self) Lipid instance.
                   point) 1 x 3 NP Array. Contains the x, y, and z 
                   coordinates, respectively, of the position 
                   axis_member) String. The name of the axis member to 
                   rotate, "Head" or "Tail".
        Requires:  extrema, xyz
        Modifies:  xyz
        Returns:   Nothing
        '''
        # Perform a rotation about each axis to align with the vector.
        # The rotation is performed by projecting along the rotation
        # axis and getting the angle between the two projections with
        # respect to the axis normal.
        # Get axis oxygen index.

        # Select which axis extreme will be mobile. 
        res  = self.extrema[axis_member][0]
        atom = self.extrema[axis_member][1]

        # z-axis rotation.
        sphere_vec  = [point[0], point[1], 0]
        protein_vec = [self.residues[res].xyz[atom, 0], 
                       self.residues[res].xyz[atom, 1], 
                       0]
        angle = gs.directional_angle(protein_vec, sphere_vec, 
                                     [0, 0, 1])
        self.rotate(axis = 'z', theta = angle)

        # y-axis rotation.
        sphere_vec  = [point[0], 0, point[2]]
        protein_vec = [self.residues[res].xyz[atom, 0], 
                       0, 
                       self.residues[res].xyz[atom, 2]]
                       
        angle = gs.directional_angle(protein_vec, sphere_vec, 
                                     [0, 1, 0])
        self.rotate(axis = 'y', theta = angle)

        # x-axis rotation.
        sphere_vec  = [0, point[1], point[2]]
        protein_vec = [0, 
                       self.residues[res].xyz[atom, 1], 
                       self.residues[res].xyz[atom, 2]]
        angle = gs.directional_angle(protein_vec, sphere_vec, 
                                     [1, 0, 0])
        self.rotate(axis = 'x', theta = angle)

    def to_pdb_file(self, atom_start = 1, resi_start = 1):
        '''
        Purpose:   Converts a Protein to a List of Strings formatted
                   as the lines of a pdb file.
        Arguments: self) Protein instance.
                   atom_start) Integer. Starting number for numbering 
                   atoms in the .pdb file.
                   resi_start) Integer. Starting number for numbering 
                   residues in the .pdb file.
        Requires:  residues
        Modifies:  Nothing
        Returns:   List of Strings. Contains the atoms of the molecule 
                   formatted as pdb lines.
        '''
        begin_atom = atom_start
        begin_resi = resi_start
        protein_lines = []
        for res in self.residues:
            protein_lines += res.to_pdb_lines(atom_start = begin_atom, 
                                              res_num    = begin_resi)
            begin_resi += 1
            begin_atom += len(res.serial)
            
        return protein_lines

    def add_mass(self):
        '''
        Purpose:   Add atomic masses to Residues of a  Protein 
                   instance.
        Arguments: self) Protein instance.
        Requires:  residues
        Modifies   residues
        Returns:   Nothing
        '''
        for res in self.residues:
            res.add_mass()

    def add_radius(self):
        '''
        Purpose:   Add atomic radii to Protein instance.
        Arguments: self) Protein instance.
        Requires:  residues
        Modifies:  residues
        Returns:   Nothing
        '''
        for res in self.residues:
            res.add_radius()

    def generate_lookup(self):
        '''
        Purpose:   Generate a dictionary for looking up an
                   atom's index for each of the Residues in
                   a Protein instance.
        Arguments: self) Protein instance.
        Requires:  residues
        Modifies:  residues
        Returns:   Nothing
        '''
        for ii in range(len(self.residues)):
            res = copy.deepcopy(self.residues[ii])
            res.generate_lookup()
            self.residues[ii] = res

if __name__ == "__main__":
    import manipulate_files as manf

    pdb = manf.PdbFile(lipid_name = "POC")
    lipid = pdb.to_lipid()
    lipid.rotate(axis = 'x', theta = 2.3)
    lines = lipid.to_pdb_lines()
    manf.write_file("test.pdb", lines)
