'''------------------------------------------------------*
| Title:   MMAEVe                                        |
|                                                        |
| Author:  Gubbin Eel (Satanic Overlord of the Swamp)    |
|                                                        |
| Purpose: Build complex Biomolecular systems for        |
|          making pretty visualizations or for molecular |
|          dynamics simulations.                         |
*-------------------------------------------------------'''

# Import required modules
import numpy as np
import pandas as pd
import scipy.spatial

'''
Shapes and Math

Code related to the distribution of points on surfaces and quaternion 
operations.
'''

def fib_lattice(length, width, height, n_points):
    ''' 
    Purpose:   Uses the fibonicci lattice to distribute points on the 
               surface of a unit square whose dimensions are then
               modified.
    Arguments: length) Float. Scale of lattice along the x-axis.
               width) Float. Scale of the lattice along the 
               y-axis.
               height) Float. z-coordinate of the points.
               n_points) Integer. Number of points to distribute.
    Returns:   N x 3 np.array of Floats. Points distributed on a
               lattice.
    '''
    if n_points == 1:
        return(np.array([[0., 0., 0.]]))
        
    eps = 0.5
    au_ratio = (1 + (5 ** 0.5)) / 2
    x = (np.arange(0, n_points, 1) / au_ratio) % 1 * length
    y = (np.arange(0, n_points, 1) + eps) \
        /                                 \
        ((n_points - 1) + (2. * eps)) * width
    z = np.repeat(0., n_points) + height
    
    return(np.stack([x, y, z], axis = -1))

def fib_disc(radius, height, n_points):
    ''' 
    Purpose:   Distribute points across the surface of a disc using 
               the fibonicci lattice.
    Arguments: radius) Float. Radius of the disc on which the points
               are distributed.
               height) Float. z-coordinate of the points.
               n_points) Integer. Number of points to distribute.
    Returns:   N x 3 np.array of Floats. Points distributed on a
               disc.
    '''
    lattice = fib_lattice(1., 1., height, n_points)
    theta = 2 * np.pi * lattice[:, 0]
    r = lattice[:, 1] ** 0.5

    x = r * radius * np.cos(theta)
    y = r * radius * np.sin(theta)
    z = lattice[:, 2]

    return(np.stack([x, y, z], axis = -1))

def fib_cylinder(radius, length, height, n_points):
    ''' 
    Purpose:   Distribute points across the surface of a cylinder 
               using the fibonicci lattice.
    Arguments: radius) Float. Radius of the cylinder.
               length) Float. Length of the cylinder.
               height) Float. Translation of the coordinates along the
               z-axis.
               n_points) Integer. Number of points to distribute.
    Returns:   N x 3 np.array of Floats. Points distributed on a
               cylinder.
    '''
    # Starting lattice
    lattice = fib_lattice(1., 1., 0., n_points)

    # Lattice => Cylindrical projection
    theta = 2 * np.pi * lattice[:, 0]
    y = lattice[:, 1]
    r = np.sqrt(y * y + 1.)

    # Cylindrical projection => Spherical coordinates
    theta = theta
    phi = np.arccos(y / r)
    r = r

    # Spherical coordinates => Cartesian coordinates
    x = r * np.sin(phi) * np.cos(theta) * radius
    y = r * np.sin(phi) * np.sin(theta) * radius
    z = r * np.cos(phi) * length

    return(np.stack([x, y, z], axis = -1))

def fib_sphere(radius, height, n_points):
    ''' 
    Purpose:   Distribute points across the surface of a sphere using 
               the fibonicc lattice.
    Arguments: radius) Float. The radius of the sphere on which points 
               are distributed.
               height) Float. Translation of the coordinates along the
               z-axis.
               num_points) Integer. The number of points to distribute 
               across the sphere's surface.
    Returns:   N x 3 np.array of Floats. Points distributed on the 
               surface of a sphere.
    '''
    # Starting lattice
    lattice = fib_lattice(1., 1., 0., n_points)

    # Lattice => Spherical projection
    theta = 2 * np.pi * lattice[:, 0]
    phi = np.arccos(1. - (2. * lattice[:, 1]))

    # Spherical projection => Cartesian coordinates
    x = np.sin(phi) * np.cos(theta) * radius
    y = np.sin(phi) * np.sin(theta) * radius
    z = np.cos(phi) * radius

    return(np.stack([x, y, z], axis = -1))

def circle(radius, height, n_points):
    '''
    Purpose:   Distribute points evenly on a circle.
    Arguments: num_points) Integer. Number of points to be 
               distributed.
               radius) Float. The desired radius of the circle.
               height) Float. Height of the system.
    Returns:   N x 3 np.array of Floats. Points distributed on a 
               circle.
    '''
    theta = np.linspace(0, 2 * np.pi, n_points, endpoint = False)
    x = radius * np.cos(theta)
    y = radius * np.sin(theta)
    z = np.repeat(height, len(x))
    
    return(np.stack([x, y, z], axis = -1))

def grid(length, width, height, l_number, w_number):
    '''
    Purpose:   Distribute l_number x w_number across a grid of 
               length x width dimensions.
    Arguments: length) Float. The desired size of the x-dimension.
               width) Float. The desired size of the y-dimension.
               height) Float. z-coordinate of the grid.
               l_number) Integer. Number of number to be distributed 
               along the x-axis.
               w_number) Integer. Number of number to be distributed 
               along the y-axis.
    Returns:   N x 3 np.array of Floats. Where each column 
               corresponds to the x, y, and z coordinates, 
               respectively.
    '''
    # Shift factor to center number in their respective grid 
    # subdivisions
    x_shift = (length / l_number) / 2.
    y_shift = (width / w_number) / 2.

    # Generate x and y coordinates
    x = np.arange(0, length, length / l_number) + x_shift
    y = np.arange(0, width, width / w_number) + y_shift
    
    # Create indicies to make combinations of x and y coordinates
    x_indicies = np.array([np.arange(0, l_number, 1)])
    x_indicies = np.repeat(x_indicies, w_number, axis = 0)
    x_indicies = x_indicies.flatten().astype(int)

    y_indicies = np.arange(0, w_number, 1)
    y_indicies = np.repeat(y_indicies, l_number, axis = 0).astype(int)
    
    # Concatenate and return points
    x = x[x_indicies]
    y = y[y_indicies]
    z = np.repeat(height, len(y_indicies))
    
    return(np.stack([x, y, z], axis = -1))

def get_rot_quat(mobile, reference):
    '''
    Purpose:   Get the rotation and inverse rotation quaternions
               needed to align a given set of vectors with another 
               set of vectors.
    Arguments: mobile) N x 3 np.array of Floats. Each row is an
               individual vector. These points are considered
               mobile.
               reference) N x 3 np.array of Floats. Each row is an
               individual point. These points will be aligned to.
    Returns:   N x 3 np.array of Floats. Where each of the original 
               points has been aligned to a vector.
    '''
    # Get the rotation axis and the angle to rotate about to align
    # "mobile" with "reference".
    rot_axis = np.cross(mobile, reference)
    dot = np.sum(mobile * reference, axis = 1)
    angle = np.arccos(dot / \
              (np.linalg.norm(mobile, axis = 1)     * \
               np.linalg.norm(reference, axis = 1)
              )
            )
    rot_axis = rot_axis / \
               np.reshape(np.linalg.norm(rot_axis, axis = 1), 
                          [len(rot_axis), 1]) # Normalize
    # Create the rotation and inverse rotation quaternions.
    cos_half = np.cos(angle / 2.)
    sin_half = np.sin(angle / 2.)
    q0, q1, q2, q3 = cos_half,                  rot_axis[:, 0] * sin_half, \
                     rot_axis[:, 1] * sin_half, rot_axis[:, 2] * sin_half
    rot_quat = np.array([q0, q1, q2, q3])
    rot_quat = rot_quat.T
    inv_quat = np.array([q0, -q1, -q2, -q3])
    inv_quat = inv_quat.T

    return(rot_quat, inv_quat)

def quat_mult(quat0, quat1):
    '''
    Purpose:   Multiply two collections of quaternions.
    Arguments: quat0) N x 4 np.array of Floats. Each row is an
               individual quaternion.
               quat1) N x 4 np.array of Floats. Each row is an
               individual quaternion.
    Returns:   N x 4 np.array of Floats. Where each row 
               corresponds to a result of the individual quaternion 
               multiplications.
    '''
    w1, x1, y1, z1 = quat0[:, 0], quat0[:, 1], quat0[:, 2], quat0[:, 3]
    w2, x2, y2, z2 = quat1[:, 0], quat1[:, 1], quat1[:, 2], quat1[:, 3]
    return(
    np.column_stack([w1*w2 - x1*x2 - y1*y2 - z1*z2,
                     w1*x2 + x1*w2 + y1*z2 - z1*y2,
                     w1*y2 - x1*z2 + y1*w2 + z1*x2,
                     w1*z2 + x1*y2 - y1*x2 + z1*w2]))

'''
Molecules

Code related to the representation and manipulation of molecules.
'''
# Atomic masses of common elements. Assumes martini parameters if no
# elemental sybol is provided. Maybe change that?
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
# Atomic radii of common elements. Assumes martini parameters if no
# elemental sybol is provided. Maybe change that?
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

class Molecule(object):
    '''
    Description: 
        Represents a single molecule. The code was developed for usage
        with lipids and proteins specifically but it should
        work for DNA, carbohydrates, small ligands, etc. as well. Else,
        adapting it to do so should require little effort. Hetero atoms
        are not currently supported as they are not considered when
        a pdb is read.

    Attributes:  
        itp_name) String. Name of the file the structure came from. 
        Used for generating a gromacs topology file.
        serial) N x 1 np.array of Integers. Serial number of the atom.
        xyz) N x 3 np.array of Floats. x, y, and z coordinates of the 
        atoms, respectively.
        mass) N x 1 np.array of Floats. Atomic masses of the atoms.
        resi) N x 1 np.array of Integers. The residue numbers of the 
        atoms.
        resn) N x 1 np.array of Strings. The residue name of the 
        atoms.
        chain) N x 1 np.array of Strings. Contains the chain IDs of 
        the atoms.
        elem) N x 1 np.array of Strings. Contains the elemental 
        symbols of the atoms.
        name) N x 1 np.array of Strings. Contains the names of the 
        atoms. 
        charge) N x 1 np.array of Floats. Contains the charge on 
        each atom.
        radius) N x 1 np.array of Floats. Contains the atomic radii of 
        the atoms.
        lookup) Dictionary. Uses unique keys for the quick lookup of 
        an atom of interest. The keys are built using the atom serial 
        numbers, residue numbers, and residue names seperated by '-'.
        Used to find the user-defined head and tail of a given 
        molecule.
        extrema) Dictionary. Meant to contain the indicies of the 
        for the user-defined head and tail of the molecule.
        residue_count) Integer. The number residues in the molecule.
        residue_atom_count) M x 1 np.array of floats. The number of 
        atoms in each of the residues.
        * N = The number of atoms in the molecule.
        * M = The number of residues in the molecule.
    '''

    def __init__(self, itp_name = "", serial = None, xyz = None, 
                 mass = None, resi = None, resn = None, chain = None, 
                 elem = None, name = None, charge = None, 
                 radius = None, lookup = None, extrema = None, 
                 residue_count = 0, residue_atom_count = None):
        ''' 
        Initialize class instance.
        '''
        self.itp_name = itp_name
        if serial is None:
            self.serial = np.array([])
        else:
            self.serial  = serial
        if xyz is None:
            self.xyz = np.array([])
        else:
            self.xyz     = xyz
        if mass is None:
            self.mass = np.array([])
        else:
            self.mass    = mass
        if resi is None:
            self.resi = np.array([])
        else:
            self.resi    = resi
        if resn is None:
            self.resn = np.array([])
        else:
            self.resn    = resn
        if chain is None:
            self.chain = np.array([])
        else:
            self.chain   = chain
        if elem is None:
            self.elem = np.array([])
        else:
            self.elem    = elem
        if name is None:
            self.name = np.array([])
        else:
            self.name    = name
        if charge is None:
            self.charge = np.array([])
        else:
            self.charge  = charge
        if radius is None:
            self.radius = np.array([])
        else:
            self.radius  = radius
        if residue_atom_count is None:
            self.residue_atom_count = np.array([]) 
        else:
            self.residue_atom_count = residue_atom_count 
        if lookup is None:
            self.lookup = dict() 
        else:
            self.lookup  = lookup
        if extrema is None:
            self.extrema = dict() 
        else:
            self.extrema = extrema

        self.residue_count = residue_count 

    def add_mass(self):
        '''
        Purpose:   Add atomic masses to a Molecule instance.
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
        self.radius = np.array([atomic_radii[elem]      \
                               for elem in self.elem ])
        
    def generate_lookup(self):
        '''
        Purpose:   Generate or update the dictionary for looking up an 
                   atom's index.
        Arguments: self) Molecule instance.
        Requires:  serial, resi, resn
        Modifies:  lookup
        Returns:   Nothing
        '''
        key = np.char.add(self.serial.astype(int).astype(str), '-')
        key = np.char.add(key, self.resi.astype(int).astype(str))
        key = np.char.add(key, '-')
        key = np.char.add(key, self.resn)
        index = np.arange(0, len(self.serial), 1)
        self.lookup = dict(zip(key, index))
    
    def determine_extrema(self, head_id, tail_id):
        '''
        Purpose:   Generate or update the dictionary containing the 
                   indicies of the molecule's head and tail.
        Arguments: self) Protein instance.
                   head_id) String. Head atom's serial number, residue
                   number, and residue name, seperated by
                   dashes. For example, '80054-911-DID', where 
                   '80054' is the serial number, '911' is the 
                   residue number, and 'DID' is the residue name.
                   tail_id) String. Tail atom's serial number, residue
                   number, and residue name, seperated by
                   dashed. For example, '421-70-VAL' where '421'
                   is the atom serial number, '70' is the 
                   residue number, and 'VAL' is the residue name.
        Requires:  lookup
        Modifies:  extrema
        Returns:   Nothing
        NOTES:     Change this such that 'lookup' is no longer needed.
                   The lookup can be performed using logical selection.
                   self.xyz[(self.serial == head_serial) and 
                            (self.resi   == head_resi)   and 
                            (self.resn   == head_resn)   and]
        '''
        self.extrema["Head"] = self.lookup[head_id]
        self.extrema["Tail"] = self.lookup[tail_id]

    def get_serial(self):
        '''
        Purpose:   Get the Molecule instance's serial numbers.
        Arguments: self) Molecule instance.
        Requires:  serial
        Modifies:  Nothing
        Returns:   N x 1 np.array of Integers. The serial attribute 
                   of the Molecule instance.
        '''
        return(self.serial)

    def get_xyz(self):
        '''
        Purpose:   Get the Molecule instance's xyz coordinates.
        Arguments: self) Molecule instance.
        Requires:  xyz
        Modifies:  Nothing
        Returns:   N x 3 np.array of Floats. The xyz attribute 
                   of the Molecule instance.
        '''
        return(self.xyz)

    def get_mass(self):
        '''
        Purpose:   Get the Molecule instance's masses.
        Arguments: self) Molecule instance.
        Requires:  mass
        Modifies:  Nothing
        Returns:   N x 1 np.array of Floats. The mass attribute 
                   of the Molecule instance.
        '''
        return(self.mass)

    def get_resi(self):
        '''
        Purpose:   Get the Molecule instance's residue numbers.
        Arguments: self) Molecule instance.
        Requires:  resi
        Modifies:  Nothing
        Returns:   N x 1 np.array of Integers. The resi attribute 
                   of the Molecule instance.
        '''
        return(self.resi)

    def get_resn(self):
        '''
        Purpose:   Get the Molecule instance's residue names.
        Arguments: self) Molecule instance.
        Requires:  resn
        Modifies:  Nothing
        Returns:   N x 1 np.array of Strings. The resi attribute 
                   of the Molecule instance.
        '''
        return(self.resn)

    def get_chain(self):
        '''
        Purpose:   Get the Molecule instance's chain ids.
        Arguments: self) Molecule instance.
        Requires:  chain
        Modifies:  Nothing
        Returns:   N x 1 np.array of Strings. The chain attribute 
                   of the Molecule instance.
        '''
        return(self.chain)

    def get_elem(self):
        '''
        Purpose:   Get the Molecule instance's element names.
        Arguments: self) Molecule instance.
        Requires:  elem
        Modifies:  Nothing
        Returns:   N x 1 np.array of Strings. The elem attribute 
                   of the Molecule instance.
        '''
        return(self.elem)

    def get_name(self):
        '''
        Purpose:   Get the Molecule instance's element names.
        Arguments: self) Molecule instance.
        Requires:  name
        Modifies:  Nothing
        Returns:   N x 1 np.array of Strings. The elem attribute 
                   of the Molecule instance.
        '''
        return(self.name)

    def get_charge(self):
        '''
        Purpose:   Get the Molecule instance's charge.
        Arguments: self) Molecule instance.
        Requires:  charge
        Modifies:  Nothing
        Returns:   N x 1 np.array of Floats. The charge attribute 
                   of the Molecule instance.
        '''
        return(self.charge)

    def get_radius(self):
        '''
        Purpose:   Get the Molecule instance's radius.
        Arguments: self) Molecule instance.
        Requires:  radius
        Modifies:  Nothing
        Returns:   N x 1 np.array of Floats. The radius attribute 
                   of the Molecule instance.
        '''
        return(self.radius)

    def get_count(self):
        '''
        Purpose:   Get the count of atomic coordinates in the system.
        Arguments: self) Molecule instance.
        Requires:  xyz
        Modifies:  Nothing
        Returns:   Integer. Number of atomic coordinates in the 
                   system. 
        '''
        return(len(self.xyz))

    def get_head(self):
        '''
        Purpose:   Get the Molecule instance's head coordinate.
        Arguments: self) Molecule instance.
        Requires:  extrema, xyz
        Modifies:  Nothing
        Returns:   N x 3 np.array of Floats. The head coordinate of
                   the Molecule instance.
        '''
        # lol
        return([self.xyz[self.extrema["Head"]]] * len(self.xyz))

    def get_tail(self):
        '''
        Purpose:   Get the Molecule instance's tail coordinate.
        Arguments: self) Molecule instance.
        Requires:  extrema, xyz
        Modifies:  Nothing
        Returns:   N x 3 np.array of Floats. The tail coordinate of
                   the Molecule instance.
        '''
        # lol
        return([self.xyz[self.extrema["Tail"]]] * len(self.xyz))

    def get_residue_count(self):
        '''
        Purpose:   Get number of residues in the molecule.
        Arguments: self) Molecule instance.
        Requires:  residue_count
        Modifies:  Nothing
        Returns:   Integer. The total number of residues in the 
                   Molecule instance.
        '''
        return(self.residue_count)

    def get_residue_atom_count(self):
        '''
        Purpose:   Get number of residues in the molecule.
        Arguments: self) Molecule instance.
        Requires:  residue_atom_count
        Modifies:  Nothing
        Returns:   M x 1 np.array of Integers. The total number of 
                   atomic coordinates in each of the Molecule 
                   instances.
        '''
        return(self.residue_atom_count)

    def get_itp_name(self):
        '''
        Purpose:   Get number of residues in the molecule.
        Arguments: self) Molecule instance.
        Requires:  residue_atom_count
        Modifies:  Nothing
        Returns:   M x 1 np.array of Integers. The total number of 
                   atomic coordinates in each of the Molecule 
                   instances.
        '''
        return(self.itp_name)

    '''
    Converting the 'Molecule.get_*' functions from pyfuncs to 
    numpy.ufuncs. This is critical for performance during system 
    initialization as it vectorized iteration over an array of 
    Molecule instances.
    '''
    get_serial = np.frompyfunc(get_serial, 1, 1)
    get_xyz    = np.frompyfunc(get_xyz, 1, 1)
    get_mass   = np.frompyfunc(get_mass, 1, 1)
    get_resi   = np.frompyfunc(get_resi, 1, 1)
    get_resn   = np.frompyfunc(get_resn, 1, 1)
    get_chain  = np.frompyfunc(get_chain, 1, 1)
    get_elem   = np.frompyfunc(get_elem, 1, 1)
    get_name   = np.frompyfunc(get_name, 1, 1)
    get_charge = np.frompyfunc(get_charge, 1, 1)
    get_radius = np.frompyfunc(get_radius, 1, 1)
    get_head   = np.frompyfunc(get_head, 1, 1)
    get_tail   = np.frompyfunc(get_tail, 1, 1)
    get_residue_count = np.frompyfunc(get_residue_count, 1, 1)
    get_residue_atom_count = np.frompyfunc(get_residue_atom_count, 
                                           1, 1)
    get_count = np.frompyfunc(get_count, 1, 1)
    get_itp_name = np.frompyfunc(get_itp_name, 1, 1)

'''
PDB File Import

Code related to the import, representation, and conversion of PDB 
files.
'''
class PdbFile(object):
    '''
    Description:
        Representation of a PdbFile. Stores all the information found 
        within a specified PDB file. The pdb file can be converted to 
        an instance of the Molecule class.

    Attributes: 
        record_name) N x 1 np.array of Strings. PDB record field.
        serial) N x 1 np.array of Integers. PDB ATOM record field: 
        Atom serial number.
        atom_name) N x 1 np.array of Strings. PDB ATOM record field: 
        Atom name.
        alternate_loc) N x 1 np.array of Strings. PDB ATOM record 
        field: Alternate location indicator.
        residue_name) N x 1 np.array of Strings. PDB ATOM record 
        field: Residue name.
        chain) N x 1 np.array of Strings. PDB ATOM record field: Chain 
        identifier.
        residue_num) N x 1 np.array of Integers. PDB ATOM record 
        field: Residue sequence number.
        insertion) N x 1 np.array of Strings. PDB ATOM record field: 
        Insertion code.
        x) N x 1 np.array of Floats. PDB ATOM record field: 
        x-coordinate.
        y) N x 1 np.array of Floats. PDB ATOM record field: 
        y-coordinate.
        z) N x 1 np.array of Floats. PDB ATOM record field: 
        z-coordinate.
        occupancy) N x 1 np.array of Floats. PDB ATOM record field: 
        Occupancy.
        temp_fac) N x 1 np.array of Floats. PDB ATOM record field: 
        Temperature factor.
        element) N x 1 np.array of Strings. PDB ATOM record field: 
        Element symbol.
        charge) N x 1 np.array of floats. PDB ATOM record field: 
        Charge.
    '''

    def __init__(self, file_name = None, structure_name = None):
        '''
        Initialize the class instance.
        '''
        ''' 
        Purpose:   Initialize an instance of the PdbFile class for
                   use. It's bad practice but it currently requires
                   that one of the optional arguments be provided 
                   else it will return an error. It will read the
                   provided file, assuming it is a PDB file and 
                   collect any required information.
        Arguments: self) PdbFile instance.
                   file_name) String. Name of the file to import with 
                   the full path, if needed (not in same directory).
                   lipid_name) Name of the lipid to import. A
                   protein_name)
        Returns:   PdbFile instance.
        '''
        # Ensure a name was provided and prep either the file_name or
        if (file_name == None) and (structure_name == None):
            error_message = "Both provided files names were of "  + \
                            "'NoneType'. Please provide a vaild " + \
                            "'file_name' or 'structure_name' to " + \
                            "initialize a 'PdbFile' instance."
            raise TypeError(error_message)
        elif (file_name != None) and (structure_name != None):
            error_message = "Provide only 'file_name' or " + \
                            "'structure_name' not both."
            raise ValueError(error_message)
        elif file_name is not None:
            itp_name = file_name.split('/')[-1]
            self.itp_name = itp_name.split('.')[0]
        elif structure_name is not None:
            self.itp_name = structure_name
            file_name = "structures/" + structure_name + ".pdb"

        # Read pdb into a DataFrame using fixed-width-file reader of 
        # pandas. Faster than having to use a Python loop to iterate 
        # through and split/format lines.
        pdb_widths = [4, 7, 5, 1,  4, 1, 4, 1, 11, 
                      8, 8, 6, 6, 10, 2, 2]
        field_types = { 0 :   str,  1 :   int,  2 :   str,  3 :   str, 
                        4 :   str,  5 :   str,  6 :   int,  7 :   str, 
                        8 : float,  9 : float, 10 : float, 11 : float, 
                       12 : float, 13 :   str, 14 :   str, 15 :   str}
        pdb_file = pd.read_fwf(file_name, widths = pdb_widths, 
                               header = None)
        pdb_file = pdb_file[pdb_file[0] == "ATOM"]
        pdb_file = pdb_file.astype(field_types)
        pdb_file.replace('nan', '', inplace = True)

        # Maybe not the cleanest way to get everything into NumPy 
        # arrays. Perhaps revisit.
        self.record_name   = np.array(pdb_file[0]).astype(str)
        self.serial        = np.array(pdb_file[1]).astype(int)
        self.atom_name     = np.array(pdb_file[2]).astype(str)
        self.alternate_loc = np.array(pdb_file[3]).astype(str)
        self.residue_name  = np.array(pdb_file[4]).astype(str)
        self.chain         = np.array(pdb_file[5]).astype(str)
        self.residue_num   = np.array(pdb_file[6]).astype(int)
        self.insertion     = np.array(pdb_file[7]).astype(str)
        self.x             = np.array(pdb_file[8]).astype(float)
        self.y             = np.array(pdb_file[9]).astype(float)
        self.z             = np.array(pdb_file[10]).astype(float)
        self.occupancy     = np.array(pdb_file[11]).astype(float)
        self.temp_fac      = np.array(pdb_file[12]).astype(float)
        self.element       = np.array(pdb_file[14]).astype(str)
        self.charge        = np.array(pdb_file[15]).astype(str)

    def to_molecule(self):
        ''' 
        Purpose:   Take a PdbFile class instance and convert it to a Protein
                   class instance.
        Arguments: self) PdbFile instance.
        Returns:   Molecule instance.
        '''
        # Seems suboptimal but was the most reliable solution I 
        # could come up with. Performance isn't a huge issue b/c
        # This conversion is only performed on the order of tens 
        # of times.

        # Get residue count and number of atoms in each residue.
        residue_atom_count = []
        residue_count = 0
        atom_count = 0
        end = 0
        current_res_num  = self.residue_num[0]
        previous_res_num = self.residue_num[0]
        for ii in range(len(self.serial)):
            current_res_num = self.residue_num[ii]
            if ((current_res_num != previous_res_num) or 
                  (end + 1 == len(self.serial))):
                if end + 1 == len(self.serial):
                    atom_count += 1
                residue_atom_count.append(atom_count)
                residue_count += 1
                atom_count = 0
            atom_count +=1
            previous_res_num = current_res_num
            end += 1
        # Join x, y, and z
        x = np.transpose(np.array([self.x]))
        y = np.transpose(np.array([self.y]))
        z = np.transpose(np.array([self.z]))
        xyz = np.concatenate([x, y, z], axis = 1)
        residue_atom_count = np.array(residue_atom_count)
        # Create molecule instance
        molecule = Molecule(itp_name = self.itp_name,
                            serial = self.serial, 
                            xyz    = xyz,  
                            resi   = self.residue_num, 
                            resn   = self.residue_name, 
                            chain  = self.chain, 
                            elem   = self.element, 
                            name   = self.atom_name, 
                            charge = self.charge,
                            residue_count = residue_count,
                            residue_atom_count = residue_atom_count)
        return(molecule)

'''
Biomolecular Complex

Code related to building biomolecular systems.
'''

def read_comp(file_name):
    ''' 
    Purpose:   Read a file and convert it to a composition dictionary.
               File should follow the following formatting convention:
    
               POPC     POPS     CHOL     POP2
               0.6      0.2      0.12     0.08
               2-1-POPC 2-1-POPC 1-1-CHOL 4-1-POP2
               7-1-POPC 7-1-POPC 7-1-CHOL 11-1-POP2
    
               Where the first row is the name of the structures that 
               will be imported. The second row is the fraction of 
               positons that each structure will occupy. The third row 
               is the 'serial-residue_number-residue_name' used to 
               specify which atom will serve as the head. The fourth 
               row is the 'serial-residue_number-residue_name' used to 
               specify which atom will serve as the tail.
    Arguments: file_name) PdbFile instance.
    Returns:   Dictionary. Properly formatted to be used when 
               initializing Lattice, Grid, Disc, Cylinder, or Sphere.
    '''
    # Read file
    with open(file_name, 'r') as f:
        contents = f.readlines()
    # Get fields and info for each individual dict
    comp = []
    names     = contents[0].split()
    fractions = contents[1].split()
    heads     = contents[2].split()
    tails     = contents[3].split()
    # Create list of dicts for each Molecule
    for ii in range(len(names)):
        entry = {}
        entry["Fraction"] = float(fractions[ii])
        entry["Head"]     = heads[ii]
        entry["Tail"]     = tails[ii]
        comp.append(entry)
    # Associate each dict with its name
    comp = dict(zip(names, comp))

    return(comp)

class BiomolComplex(object):
    ''' 
    Description:
        Generic parent class used to define basic features necessary
        for or common to the building/manipulation of more specific
        biomolecular structures or components of biomolecular 
        structures.
        
    Attributes: 
        positions) (N, 3) np.array of Floats. Reference points used 
        to place individual proteins.
        structures) (M, ) np.array of Molecule instances. Used to 
        populate the attributes related to the protein structure.
        ratios) (M, ) np.array of Floats ≤ 1 whose sum is 1.0. 
        Fraction of 'positions' that should correspond to each member 
        of 'structures'.
        ids) (N, ) np.array of Integers. Each integer represents an 
        index corresponding to the member of structures that 
        should be placed at the corresponding member of positions.
        serial) (D, ) np.array of Integers. Atomic serial numbers 
        collected from Molecule instances.
        name) (D, ) np.array of Strings. Atom names collected from 
        Molecule instances.
        resi) (D, ) np.array of Integers. Residue numbers collected 
        from Molecule instances.
        resn) (D, ) np.array of Strings. Residue names collected from 
        Molecule instances.
        chain) (D, ) np.array of Strings/Integers. Chain IDs collected 
        from Molecule instances.
        xyz) (D, 3) np.array of Floats. Cartesian coordinates 
        collected from Molecule instances.
        elem) (D, ) np.array of Strings. Elemental symbols collected 
        from  instances.
        head) (D, 3) np.array of Floats. An array with coordinates 
        of the defined Molecule head.
        tail) (D, 3) np.array of Floats. An array with coordinates 
        of the defined Molecule tail.
        entity_count) Integer. Number of total entities. 
        entity_atom_count) (N, ) np.array of Integers. Number of atoms 
        in each entity. 
        residue_count) (N, ) np.array of Integers. Number of residues 
        in each entity.
        residue_atom_count) (E, ) np.array of Integers. Number of atoms 
        in each residue.
        * N) Number of positions Molecules will be assigned to.
        * M) Number of Molecule types they system is composed of.
        * D) Number of atomic coordinates in the system.
        * E) Number of residues in the system.
    '''
    def __init__(self, positions = None, structures = None, 
                 ratios = None, ids = None, serial = None, 
                 name = None, resi = None, resn = None, chain = None, 
                 xyz = None, elem = None, mass = None, charge = None, 
                 radius = None, head = None, tail = None, 
                 entity_count = 0, entity_atom_count = None, 
                 residue_count = None, residue_atom_count = None):
        ''' 
        Initialize class instance.
        '''
        if positions is None:
            self.positions  = np.array([]) 
        else:
            self.positions = positions 
        if structures is None:
            self.structures = np.array([]) 
        else:
            self.structures = structures 
        if ratios is None:
            self.ratios = np.array([]) 
        else:
            self.ratios = ratios 
        if ids is None:
            self.ids = np.array([])
        else:
            self.ids = ids 
        if serial is None: 
            self.serial = np.array([])
        else:
            self.serial = serial 
        if name is None: 
            self.name = np.array([])
        else:
            self.name = name 
        if resi is None: 
            self.resi = np.array([])
        else:
            self.resi = resi 
        if resn is None: 
            self.resn = np.array([])
        else:
            self.resn = resn 
        if chain is None: 
            self.chain = np.array([])
        else:
            self.chain = chain 
        if xyz is None: 
            self.xyz = np.array([])
        else:
            self.xyz = xyz 
        if elem is None: 
            self.elem = np.array([])
        else:
            self.elem = elem 
        if mass is None: 
            self.mass = np.array([])
        else:
            self.mass = mass 
        if charge is None: 
            self.charge = np.array([])
        else:
            self.charge = charge 
        if radius is None: 
            self.radius = np.array([])
        else:
            self.radius = radius 
        if head is None: 
            self.head = np.array([])
        else:
            self.head = head 
        if tail is None: 
            self.tail = np.array([])
        else:
            self.tail = tail 
        if entity_atom_count is None: 
            self.entity_atom_count = np.array([])
        else:
            self.entity_atom_count = entity_atom_count 
        if residue_count is None: 
            self.residue_count = np.array([])
        else:
            self.residue_count = residue_count 
        if residue_atom_count is None: 
            self.residue_atom_count = np.array([])
        else:
            self.residue_atom_count = residue_atom_count 
        self.entity_count  = entity_count

    def __add__(self, struc2):
        '''
        Purpose:   Addition operator override. Used to translate
                   a system by a fixed 3D vector or concatenate two
                   systems.
        Arguments: self) lipidStructure instance.
                   struc2) (3, ) np.array of Floats.
                              --- OR ---
                   struc2) BiomolComplex instance.
        Requires:  xyz, head, tail, positions
                   --- OR ---
                   all
        Modifies:  xyz, head, tail, positions
                   --- OR ---
                   Nothing
        Returns:   Nothing
                   --- OR ---
                   BiomolComplex instance.
        * Probably asking a little too much of the humble addition 
          operator but having experimented with others for 
          concatenation I can say that none of them _feel_ correct. 
          And other operators for simple translation would just be 
          sacrilige. Probably need to revisit.
        '''
        if (isinstance(struc2, np.ndarray)) and (struc2.shape == (3,)):
            self.xyz += struc2
            self.head += struc2
            self.tail += struc2
            self.positions += struc2
        elif isinstance(struc2, BiomolComplex):
            return(self.concatenate(struc2))
        else:
            error_message = "Addition called for an unsupported " + \
                            "type. Please,ensure that the "       + \
                            "operation is called with a "         + \
                            "suppored type."
            raise ValueError(error_message)

    def concatenate(self, struc2):
        '''
        Purpose:   Concatenate fields of two instance of BiomolComplex 
                   to create a new BiomolComplex.
        Arguments: self) BimolComplex instance.
                   struc2) BimolComplex instance.
        Requires:  all
        Modifies:  Nothing
        Returns:   BiomolComplex 
        * Little messy. Maybe add a loop or create a vectorized 
          equivalent? Is it really worth it?
        '''
        # Bunch of np.concatenate calls
        positions          = np.concatenate([self.positions, 
                                             struc2.positions])
        structures         = np.concatenate([self.structures, 
                                             struc2.structures])
        serial             = np.concatenate([self.serial, 
                                             struc2.serial])
        name               = np.concatenate([self.name, 
                                             struc2.name])
        resi               = np.concatenate([self.resi, 
                                             struc2.resi])
        resn               = np.concatenate([self.resn, 
                                             struc2.resn])
        chain              = np.concatenate([self.chain, 
                                             struc2.chain])
        xyz                = np.concatenate([self.xyz, 
                                             struc2.xyz])
        elem               = np.concatenate([self.elem, 
                                             struc2.elem])
        head               = np.concatenate([self.head, 
                                             struc2.head])
        tail               = np.concatenate([self.tail, 
                                             struc2.tail])
        entity_atom_count  = np.concatenate([self.entity_atom_count, 
                                             struc2.entity_atom_count])
        residue_count      = np.concatenate([self.residue_count, 
                                             struc2.residue_count])
        residue_atom_count = np.concatenate([self.residue_atom_count, 
                                             struc2.residue_atom_count])
        mass               = np.concatenate([self.mass, 
                                             struc2.mass])
        charge             = np.concatenate([self.charge, 
                                             struc2.charge])
        radius             = np.concatenate([self.radius, 
                                             struc2.radius])

        # Modify and combine ratios and entity_counts
        ratios = np.concatenate([self.ratios * self.entity_count, 
                                 struc2.ratios * struc2.entity_count])
        entity_count = self.entity_count + struc2.entity_count
        ratios = ratios / entity_count

        # Update the ids to reflect the new, combined structure list.
        ids = np.concatenate([self.ids, 
                              struc2.ids + len(self.structures)])

        return(
            BiomolComplex(positions = positions, 
                structures = structures, ratios = ratios, ids = ids, 
                serial = serial, name = name, resi = resi, 
                resn = resn, chain = chain, xyz = xyz, elem = elem, 
                mass = mass, charge = charge, radius = radius, 
                head = head, tail = tail, entity_count = entity_count, 
                entity_atom_count = entity_atom_count, 
                residue_count = residue_count, 
                residue_atom_count = residue_atom_count))

    def collect_blocks(self, composition = dict()):
        '''
        Purpose:   To import and store specified files as Molecule 
                   instances and get their corresponding ratios.
        Arguments: self) BiomolComplex instance.
                   composition) Dictionary. Key is String and Value is 
                   another dictionary. The nested dictionary should 
                   have three Keys "Fraction", "Head", "Tail". The 
                   Keys of the main ditionary are the names of the 
                   proteins that you want to use. For the secondary 
                   Dictionary, "Fraction" should be a Float <= 1.0 
                   with all "Fractions" summing to 1.0. "Head" and 
                   "Tail" are strings. They should be formatted as 
                   described in the Molecule.determine_extrema method 
                   An example is provided below.
                   composition = {"APO" : {"Fraction" : 0.5 
                                           "Head"     : "123-12-ARG"
                                           "Tail"     : "444-87-GLY"}
                                  "ANX" : {"Fraction" : 0.5
                                           "Head"     : "336-38-TYR"
                                           "Tail"     : "901-99-ASP"}
                                 }
        Requires:  Nothing
        Modifies:  structures, ratios
        Returns:   Nothing
        * See the 'get_composition_dict' for a function to get the
          composition from a file.
        '''
        self.ratios = []
        self.structures = []
        for name in composition:
            self.ratios.append(composition[name]["Fraction"])
            pdb = PdbFile(structure_name = name)
            mol = pdb.to_molecule()
            mol.generate_lookup()
            mol.determine_extrema(composition[name]["Head"], 
                                  composition[name]["Tail"])
            mol.add_radius()
            mol.add_mass()
            self.structures.append(mol)
        self.structures = np.array(self.structures)
        self.ratios = np.array(self.ratios)

    def set_seed(seed):
        np.random.seed(seed)

    def gen_ids(self, seed = None):
        ''' 
        Purpose:   Assign an ID for each point generated by the 
                   'generate_positions' method. 'generate_positions'
                   is defined for the child classes of BiomolComplex.
                   These IDs correspond to which Molecule instance 
                   will be placed at the corresponding point.
        Arguments: self) BiomolComplex instance.
                   seed) Integer. Random seed that can be used to 
                   produce consisten results across runs.
        Requires:  structures, positions, ratios
        Modifies:  ids
        Returns:   Nothing
        '''
        if seed == None:
            state = np.random.RandomState()
        else:
            state = np.random.RandomState(seed)
        id_pool = np.arange(0, len(self.structures), 1)
        self.ids = state.choice(id_pool, len(self.positions),
                                p = self.ratios)

    def sort_ids(self):
        ''' 
        Purpose:   Sort the ids with their corresponding positions 
                   such that all instances occur in sequence. This 
                   is necessary for generating a gromacs topology 
                   with a reasonable number of entries.
        Arguments: self) BiomolComplex instance.
        Requires:  ids, positions
        Modifies:  positions
        Returns:   Nothing
        '''
        sort_index = np.argsort(self.ids)
        self.ids = self.ids[sort_index]
        self.positions = self.positions[sort_index]

    def populate_structure_components(self):
        ''' 
        Purpose:   Get starting values for all properties of the 
                   system.
        Arguments: self) BiomolComplex instance.
        Requires:  ids, structures
        Modifies:  serial, xyz, mass, resi, resn, chain, elem, name, 
                   charge, radius, head, tail, residue_count, 
                   residue_atom_count, entity_count, entity_atom_count
        Returns:   Nothing
        '''
        # Generate np.array of structure instances
        molecules = self.structures[self.ids]

        # Mostly calls to the vectorize "Molecule.get_*" methods
        self.serial = np.concatenate(Molecule.get_serial(molecules))
        self.serial = self.serial.astype(int)
        self.xyz    = np.concatenate(Molecule.get_xyz(molecules))
        self.mass   = np.concatenate(Molecule.get_mass(molecules))
        self.resi   = np.concatenate(Molecule.get_resi(molecules))
        self.resi   = self.resi.astype(int)
        self.resn   = np.concatenate(Molecule.get_resn(molecules))
        self.chain  = np.concatenate(Molecule.get_chain(molecules))
        self.elem   = np.concatenate(Molecule.get_elem(molecules))
        self.name   = np.concatenate(Molecule.get_name(molecules))
        self.charge = np.concatenate(Molecule.get_charge(molecules))
        self.radius = np.concatenate(Molecule.get_radius(molecules))
        self.head   = np.concatenate(Molecule.get_head(molecules))
        self.tail   = np.concatenate(Molecule.get_tail(molecules))

        self.residue_count = Molecule.get_residue_count(molecules)
        self.residue_count = self.residue_count.astype(int)

        self.residue_atom_count = \
            np.concatenate(Molecule.get_residue_atom_count(molecules))

        self.residue_atom_count = self.residue_atom_count.astype(int)
        self.entity_count       = len(self.positions)
        self.entity_atom_count  = Molecule.get_count(molecules)
        self.entity_atom_count  = self.entity_atom_count.astype(int)

    def regenerate_residue_info(self):
        ''' 
        Purpose:   Update resiude information after removing 
                   entities from the system as is done in the
                   BiomolComplex.remove_overlap method.
        Arguments: self) BiomolComplex instance.
        Requires:  ids, structures
        Modifies:  residue_count, residue_atom_count
        Returns:   Nothing
        '''
        molecules = self.structures[self.ids]

        self.residue_count      = Molecule.get_residue_count(molecules)
        self.residue_count      = self.residue_count.astype(int)
        self.residue_atom_count = \
            np.concatenate(Molecule.get_residue_atom_count(molecules))
        self.residue_atom_count = self.residue_atom_count.astype(int)

    def regenerate_entity_info(self):
        ''' 
        Purpose:   Update entity information after removing 
                   entities from the system as is done in the
                   BiomolComplex.remove_overlap method.
        Arguments: self) BiomolComplex instance.
        Requires:  ids, structures, positions
        Modifies:  entity_count, entity_atom_count
        Returns:   Nothing
        '''
        molecules = self.structures[self.ids]

        self.entity_count      = len(self.positions)
        self.entity_atom_count = Molecule.get_count(molecules)
        self.entity_atom_count = self.entity_atom_count.astype(int)

    def get_entity_centroids(self):
        def get_centroid(a, b):
            return(np.mean(self.xyz[a:b], axis = 0))

        get_centroid = np.frompyfunc(get_centroid, 2, 1)

        entity_sele = np.cumsum(self.entity_atom_count)
        entity_sele = np.concatenate([np.array([0]), 
                                       entity_sele])
        lower_sele = entity_sele[:-1]
        upper_sele = entity_sele[1:]

        return(np.vstack(get_centroid(lower_sele, upper_sele)))


    def remove_overlap(self, struc2, radius, atom_based = True):
        ''' 
        Purpose:   Remove entities in self that overlap with entities 
                   in struct2.
        Arguments: self) BiomolComplex instance.
                   struc2) BiomolComplex instance.
                   radius) Float. Distance cuttoff for atoms of self 
                   and struc2 which will trigger removal of self 
                   entity.        
                   atom_based) Boolean. Controls whether atomic or 
                   entity centroids are used for removing overlap.
        Requires:  ids, structures, positions
        Modifies:  all
        Returns:   Nothing
        '''

        struc1_tree = scipy.spatial.KDTree(self.xyz)
        if atom_based:
            struc2_tree = scipy.spatial.KDTree(struc2.xyz)
        else:
            struc2_tree = \
                scipy.spatial.KDTree(struc2.get_entity_centroids())

        sparse = \
            struc1_tree.sparse_distance_matrix(struc2_tree, radius, 
                                               output_type = "ndarray")

        # Remove lipids that overlap with proteins
        hit_index = list(set([dist[0] for dist in sparse]))

        # Code works well but is super nasty. Clean up and add 
        # comments.
        if (len(hit_index) > 0):
            removal_indicies = \
                np.repeat(np.arange(0, self.entity_count, 1), 
                                    self.entity_atom_count)
            removal_indicies = list(set(removal_indicies[hit_index]))

            removal_sele = np.cumsum(self.entity_atom_count)
            removal_sele = np.concatenate([np.array([0]), 
                                           removal_sele])
            lower_sele = removal_sele[:-1]
            upper_sele = removal_sele[1:]
            lower_sele = lower_sele[removal_indicies]
            upper_sele = upper_sele[removal_indicies]
            removal_sele = np.concatenate(
                [np.arange(lower_sele[ii], upper_sele[ii], 1) \
                 for ii in range(len(lower_sele))]            \
            )

            self.serial = np.delete(self.serial, removal_sele)
            self.xyz    = np.delete(self.xyz, removal_sele, axis = 0)
            self.mass   = np.delete(self.mass, removal_sele)
            self.resi   = np.delete(self.resi, removal_sele)
            self.resn   = np.delete(self.resn, removal_sele)
            self.chain  = np.delete(self.chain, removal_sele)
            self.elem   = np.delete(self.elem, removal_sele)
            self.name   = np.delete(self.name, removal_sele)
            self.head   = np.delete(self.head, removal_sele, axis = 0)
            self.tail   = np.delete(self.tail, removal_sele, axis = 0)
            self.charge = np.delete(self.charge, removal_sele)
            self.radius = np.delete(self.radius, removal_sele)

            self.ids        = np.delete(self.ids, removal_indicies)
            self.positions  = np.delete(self.positions, 
                                        removal_indicies, axis = 0)

            struc_sele      = np.array(list(set(self.ids)))
            self.structures = self.structures[struc_sele]

            self.regenerate_residue_info()
            self.regenerate_entity_info()

    def to_origin(self):
        '''
        Purpose:   Translates all of the atoms in the system to the
                   origin using the tail atom.
        Arguments: self) BiomolComplex instance.
        Requires:  head 
        Modifies:  xyz, head, tail
        Returns:   Nothing
        '''
        # Translate tail to origin
        translate = [0., 0., 0.] - self.tail
        self.xyz  = self.xyz + translate
        self.head = self.head + translate
        self.tail = self.tail + translate

    def rot_quat_mult(self, rot_quat, inv_quat):
        '''
        Purpose:   Rotates the system according to defined a pre-
                   defined rotation and inverse rotation quaternion.
        Arguments: self) BiomolComplex instance.
                   rot_quat) (D, 4) np.array of Floats. Rotation 
                   quaternion for each indiviudal atom of the system.
                   inv_quat) (D, 4) np.array of Floats. Inverse 
                   rotation quaternion for each indiviudal atom of the 
                   system.
        Requires:  xyz, head, tail 
        Modifies:  xyz, head, tail
        Returns:   Nothing
        '''
        # Transform xyz to quaternions
        self.xyz  = np.insert(self.xyz, 0, 0., axis = 1)
        self.head = np.insert(self.head, 0, 0., axis = 1)
        self.tail = np.insert(self.tail, 0, 0., axis = 1)

        # Perform quaternion multiplication
        self.xyz = quat_mult(quat_mult(rot_quat, self.xyz), 
                             inv_quat)[:, 1:]
        self.head = quat_mult(quat_mult(rot_quat, self.head), 
                              inv_quat)[:, 1:]
        self.tail = quat_mult(quat_mult(rot_quat, self.tail), 
                              inv_quat)[:, 1:]

    def to_positions(self, positions):
        '''
        Purpose:   Translates all of the atoms in the system to their
                   positions using the head atom.
        Arguments: self) BiomolComplex instance.
                   positions) (D, ) np.array of Floats. The positions 
                   that each individual atom should be translated to.
        Requires:  head 
        Modifies:  xyz, head, tail
        Returns:   Nothing
        '''
        # Translate head to position
        translate  = positions - self.head
        self.xyz  += translate
        self.head += translate
        self.tail += translate

    def centroid(self):
        '''
        Purpose:   Calculate the centroid of a system.
        Arguments: self) BiomolComplex instance.
        Requires:  xyz
        Modifies:  none
        Returns:   (3, ) np.array of Floats. System centroid.
        '''
        return(np.sum(self.xyz, axis = 0) / len(self.xyz))

    def min(self):
        '''
        Purpose:   Determine the minimum xyz coords of the system.
        Arguments: self) BiomolComplex instance.
        Requires:  xyz
        Modifies:  none
        Returns:   (3, ) np.array of Floats. Min xyz coords.
        '''
        return(np.min(self.xyz, axis = 0))

    def max(self):
        '''
        Purpose:   Determine the maximum xyz coords of the system.
        Arguments: self) BiomolComplex instance.
        Requires:  xyz
        Modifies:  none
        Returns:   (3, ) np.array of Floats. Max xyz coords.
        '''
        return(np.max(self.xyz, axis = 0))

    def dim(self):
        '''
        Purpose:   Determine the dimensions of the system.
        Arguments: self) BiomolComplex instance.
        Requires:  xyz
        Modifies:  none
        Returns:   (3, ) np.array of Floats. System dimensions.
        '''
        return(np.abs(self.max - self.min))


    def swap_yz(self):
        '''
        Purpose:   Swap the y and z coordinates. This is nice because 
                   different software adhere to different conventions 
                   for which axis is considered 'height'.
        Arguments: self) BiomolComplex instance.
                   file_name) String. Name to write the top file to.
                   should include the '.top' extension in the name.
        Requires:  xyz, head, tail, positions
        Modifies:  xyz, head, tail, positions
        Returns:   Nothing
        '''
        tmp = self.xyz[:, 1]
        self.xyz[:, 1] = self.xyz[:, 2]
        self.xyz[:, 2] = tmp

        tmp = self.head[:, 1]
        self.head[:, 1] = self.head[:, 2]
        self.head[:, 2] = tmp

        tmp = self.tail[:, 1]
        self.tail[:, 1] = self.tail[:, 2]
        self.tail[:, 2] = tmp

        tmp = self.positions[:, 1]
        self.positions[:, 1] = self.positions[:, 2]
        self.positions[:, 2] = tmp

    def write_cif(self, file_name, keep_asym = True, buff = 8192, 
                  chunksize = 100000):
        '''
        Purpose:   Write structure to a cif file. No advanced field 
                   manipulation is supported. This is a bare-bones 
                   method to get the structure written as a cif file. 
                   I have had spotty success with getting the 
                   structures read. PyMol will read them without 
                   issue. Meant for a quick visualization.
        Arguments: self) BiomolComplex instance.
                   file_name) String. Name to write the cif file to.
                   should include the '.cif' extension in the name.
                   keep_asym) Boolean. Controls whether an attempt 
                   to keep the original asym id is made or if a new 
                   one is generated.
                   buff) Integer. Controls the buffersize of the open
                   function. 
                   chunksize) Integer. Controls the number of "chunks" 
                   to write the file in. 
        Requires:  residue_count, tail, name, residue_atom_count, 
                   resn, xyz, name, chain 
        Modifies:  none
        Returns:   Nothing
        '''
        atom_count = len(self.tail)            # Total number of atoms
        res_count = np.sum(self.residue_count) # Total # of residues

        if keep_asym == True:
            # Keeps original, replaces any blank fields
            label_asym_id = self.chain
            label_asym_id[label_asym_id == ''] = '+'
            label_asym_id = np.char.replace(self.chain, '+', '?')
            pass
        else:
            # Makes asym_id uniform
            label_asym_id = np.repeat(np.array(["A"]), atom_count)

        # PDB fields to write: atom serial, atom name, residue name, 
        # chain id, residue number, x, y, z.
        fields = [
            np.mod(np.arange(1, atom_count + 1, 1), 100000), # len correctable
            self.name,
            self.resn,
            label_asym_id,
            np.mod(np.repeat(np.arange(1, res_count + 1, 1), 
                          self.residue_atom_count), 10000), # len correctable
            np.round(self.xyz[:, 0], decimals = 3),
            np.round(self.xyz[:, 1], decimals = 3),
            np.round(self.xyz[:, 2], decimals = 3)]
        nparts = self.xyz.shape[0] # Number of file lines
        resid = np.zeros(nparts, dtype = int)
        count = 1

        # Used to format when writing file lines.
        write_string = "ATOM {:5d} {:4s} {:4s} {:1s} {:4d} " + \
                       "{:8.3f} {:8.3f} {:8.3f}\n"
        with open(file_name, 'w', buff) as fout:
            chunks = list(range(0, nparts, chunksize))
            chunks.append(nparts + 1)   # so slicing is from 2ndtolast:size
            header = "data_LIPID\n#\nloop_\n_atom_site.group_PDB\n" + \
                     "_atom_site.id\n_atom_site.label_atom_id\n"    + \
                     "_atom_site.label_comp_id\n"                   + \
                     "_atom_site.label_asym_id\n"                   + \
                     "_atom_site.label_seq_id\n"                    + \
                     "_atom_site.Cartn_x\n_atom_site.Cartn_y\n"     + \
                     "_atom_site.Cartn_z\n"
            fout.writelines(header)

            for startind, stopind in zip(chunks[:-1], chunks[1:]):
                fout.writelines(
                    [write_string.format(ii[0], ii[1], ii[2], ii[3], 
                                         ii[4], ii[5], ii[6], ii[7]) 
                     for ii in zip(fields[0][startind:stopind],  
                                   fields[1][startind:stopind],
                                   fields[2][startind:stopind],
                                   fields[3][startind:stopind],
                                   fields[4][startind:stopind],
                                   fields[5][startind:stopind],
                                   fields[6][startind:stopind],
                                   fields[7][startind:stopind])
                    ]
                )

    def write_pdb(self, file_name, buff = 8192, chunksize = 100000):
        '''
        Purpose:   Write structure to a pdb file.
        Arguments: self) BiomolComplex instance.
                   file_name) String. Name to write the pdb file to.
                   should include the '.pdb' extension in the name.
                   buff) Integer. Controls the buffersize of the open
                   function. 
                   chunksize) Integer. Controls the number of "chunks" 
                   to write the file in. 
        Requires:  residue_count, tail, name, resn, chain, 
                   residue_atom_count, xyz, entity_atom_count
        Modifies:  none
        Returns:   Nothing
        '''
        atom_count = len(self.tail)            # Total number of atoms
        res_count = np.sum(self.residue_count) # Total # of residues
        write_string = "{:4s}  {:5d} {:4s} {:4s}{:1s}{:4d}    " + \
                       "{:8.3f}{:8.3f}{:8.3f}\n"

        # PDB fields to write: atom serial, atom name, residue name, 
        # chain id, residue number, x, y, z.
        fields = [
            np.repeat("ATOM", atom_count),
            np.arange(1, atom_count + 1, 1), # len correctable
            self.name,
            self.resn,
            self.chain,
            np.repeat(np.arange(1, res_count + 1, 1), 
                          self.residue_atom_count), # len correctable
            np.round(self.xyz[:, 0], decimals = 3),
            np.round(self.xyz[:, 1], decimals = 3),
            np.round(self.xyz[:, 2], decimals = 3),
            np.repeat(write_string, atom_count)]
        # Reduce atom serial and residue number to proper size
        fields[1] = np.mod(fields[1], 100000)
        fields[5] = np.mod(fields[5], 10000)

        # Join fields together into single array and insert TER lines
        insertion_indicies = np.cumsum(self.entity_atom_count)
        fields[0] = np.insert(fields[0], insertion_indicies, "TER", 
                              axis = 0)
        fields[1] = np.insert(fields[1], insertion_indicies, 0, axis = 0)
        fields[2] = np.insert(fields[2], insertion_indicies, '', axis = 0)
        fields[3] = np.insert(fields[3], insertion_indicies, '', axis = 0)
        fields[4] = np.insert(fields[4], insertion_indicies, '', axis = 0)
        fields[5] = np.insert(fields[5], insertion_indicies, 0, axis = 0)
        fields[6] = np.insert(fields[6], insertion_indicies, 0., axis = 0)
        fields[7] = np.insert(fields[7], insertion_indicies, 0., axis = 0)
        fields[8] = np.insert(fields[8], insertion_indicies, 0., axis = 0)
        write_ter = "{:3s}\n"
        fields[9] = np.insert(fields[9], insertion_indicies, write_ter, axis = 0)

        # Used to format when writing file lines.
        nparts = len(fields[0]) # Number of file lines
        count = 1
        with open(file_name, 'w', buff) as fout:
            chunks = list(range(0, nparts, chunksize))
            chunks.append(nparts + 1)   # so slicing is from 2ndtolast:size

            for startind, stopind in zip(chunks[:-1], chunks[1:]):
                fout.writelines(
                    [ii[9][startind:stopind].format(ii[0], ii[1], ii[2], ii[3], 
                                         ii[4], ii[5], ii[6], ii[7],
                                         ii[8]) 
                     for ii in zip(fields[0][startind:stopind],  
                                   fields[1][startind:stopind],
                                   fields[2][startind:stopind],
                                   fields[3][startind:stopind],
                                   fields[4][startind:stopind],
                                   fields[5][startind:stopind],
                                   fields[6][startind:stopind],
                                   fields[7][startind:stopind],
                                   fields[8][startind:stopind],
                                   fields[9][startind:stopind])
                    ]
                )

    def write_gromacs_top(self, file_name):
        '''
        Purpose:   Write a gromacs '.top' file for the structure.
        Arguments: self) BiomolComplex instance.
                   file_name) String. Name to write the top file to.
                   should include the '.top' extension in the name.
        Requires:  structures
        Modifies:  none
        Returns:   Nothing
        '''
        itp_names = Molecule.get_itp_name(self.structures)
        itp_names = itp_names.astype(str)
        _, counts = np.unique(self.ids, return_counts = True)
        counts = counts.astype(str)
        top_lines = np.char.add(itp_names, '\t\t')
        top_lines = np.char.add(top_lines, counts)
        top_header = ["; You'll need to add all of the itp-file " + \
                      "references", " ", "[ molecule ]", 
                      ";Name\t\tCount"]
        top_lines = np.concatenate([top_header, top_lines])

        np.savetxt(file_name, top_lines, fmt = "%s")

class Lattice(BiomolComplex):
    def __init__(self, length, width, height, number, 
                 composition, seed = None):
        super().__init__(self)
        self.collect_blocks(composition = composition)
        self.generate_positions(length, width, height, number)
        self.gen_ids(seed = seed)
        self.sort_ids()
        self.populate_structure_components()

    def generate_positions(self, length, width, height, number):
        '''
        Purpose:   Distribute points across a plane.
        Arguments: self) BiomolComplex instance.
                   length) Float. Size of the plane in the x-dim.
                   width) Float. Size of the plane in the y-dim.
                   height) Float. Amount to shift the z-direction.
                   number) Integer. Number of positions to distribute 
                   across the surface of the plane.
        Requires:  None
        Modifies:  positions
        Returns:   None
        '''
        self.positions = fib_lattice(length, width, height, number)

    def distribute(self):
        '''
        Purpose:   Orient and place molecules at their asigned 
                   positions. The procedure translates the molecules 
                   to the origin using the tail atom. Aligns the 
                   molecules to the z-axis. Distributes the molecules 
                   to their positions using the head atoms.
        Arguments: self) BiomolComplex instance.
        Requires:  positions, entity_atom_count, xyz, head, tail
        Modifies:  xyz, head, tail
        Returns:   None
        '''
        positions = np.repeat(self.positions, self.entity_atom_count, 
                              axis = 0)

        # Translate tail to origin
        self.to_origin()

        #  Rotate head to position
        z_axis = np.repeat(np.array([[0., 0., 1.]]), len(self.xyz),
                           axis = 0)
        rot_quat, inv_quat = get_rot_quat(self.head, z_axis)
        self.rot_quat_mult(rot_quat, inv_quat)

        # Translate head to position
        self.to_positions(positions)

class Grid(Lattice):
    def __init__(self, length, width, height, l_number, w_number,
                 composition, seed = None):
        super(BiomolComplex).__init__()
        self.collect_blocks(composition = composition)
        self.generate_positions(length, width, height, l_number, 
                                w_number)
        self.gen_ids(seed = seed)
        self.sort_ids()
        self.populate_structure_components()

    def generate_positions(self, length, width, height, l_number, 
                           w_number):
        '''
        Purpose:   Distribute points across a plane.
        Arguments: self) BiomolComplex instance.
                   length) Float. Size of the plane in the x-dim.
                   width) Float. Size of the plane in the y-dim.
                   height) Float. Amount to shift the z-direction.
                   number) Integer. Number of positions to distribute 
                   across the surface of the plane.
        Requires:  None
        Modifies:  positions
        Returns:   None
        '''
        self.positions = grid(length, width, height, l_number, 
                              w_number)

class Disc(Lattice):
    def __init__(self, radius, height, number, composition, 
                 seed = None):
        super(BiomolComplex).__init__()
        self.collect_blocks(composition = composition)
        self.generate_positions(radius, height, number) 
        self.gen_ids(seed = seed)
        self.sort_ids()
        self.populate_structure_components()

    def generate_positions(self, radius, height, number):
        '''
        Purpose:   Distribute points across a disc.
        Arguments: self) BiomolComplex instance.
                   radius) Float. Radius of the disc points are to 
                   be distributed on.
                   height) Float. Amount to shift the z-direction.
                   number) Integer. Number of positions to distribute 
                   across the surface of the disc.
        Requires:  None
        Modifies:  positions
        Returns:   None
        '''
        self.positions  = fib_disc(radius, height, number)

class Circle(Lattice):
    def __init__(self, radius, height, number, composition, 
                 seed = None):
        super(BiomolComplex).__init__()
        self.collect_blocks(composition = composition)
        self.generate_positions(radius, height, number) 
        self.gen_ids(seed = seed)
        self.sort_ids()
        self.populate_structure_components()

    def generate_positions(self, radius, height, number):
        '''
        Purpose:   Distribute points across a disc.
        Arguments: self) BiomolComplex instance.
                   radius) Float. Radius of the disc points are to 
                   be distributed on.
                   height) Float. Amount to shift the z-direction.
                   number) Integer. Number of positions to distribute 
                   across the surface of the disc.
        Requires:  None
        Modifies:  positions
        Returns:   None
        '''
        self.positions  = circle(radius, height, number)

class Cylinder(BiomolComplex):
    def __init__(self, radius, length, height, number, 
                 composition, seed = None):
        super().__init__(self)
        self.collect_blocks(composition = composition)
        self.generate_positions(radius, length, height, number)
        self.gen_ids(seed = seed)
        self.sort_ids()
        self.populate_structure_components()

    def generate_positions(self, radius, length, height, number):
        '''
        Purpose:   Distribute points across a cylinder.
        Arguments: self) BiomolComplex instance.
                   radius) Float. Radius of the cylinder.
                   length) Float. How 'tall' the cylinder is.
                   height) Float. Amount to shift the z-direction.
                   number) Integer. Number of positions to distribute 
                   across the surface of the cylinder.
        Requires:  None
        Modifies:  positions
        Returns:   None
        '''
        self.positions = fib_cylinder(radius, length, height, 
                                         number)

    def distribute(self):
        '''
        Purpose:   Orient and place molecules at their asigned 
                   positions. The procedure translates the molecules 
                   to the origin using the tail atom. Aligns the 
                   molecules to the xy projection of their assigned 
                   points. Distributes the molecules to their 
                   positions using the head atoms.
        Arguments: self) BiomolComplex instance.
        Requires:  positions, entity_atom_count, xyz, head, tail
        Modifies:  xyz, head, tail
        Returns:   None
        '''
        positions = np.repeat(self.positions, self.entity_atom_count, 
                              axis = 0)

        # Translate tail to origin
        self.to_origin()

        # Align to xy projection
        xy_proj = positions[:, :2]
        xy_proj = np.insert(xy_proj, 2, 0., axis = 1)
        rot_quat, inv_quat = get_rot_quat(self.head, xy_proj)
        self.rot_quat_mult(rot_quat, inv_quat)

        # Translate head to position
        self.to_positions(positions)

class Sphere(BiomolComplex):
    ''' 
    Description:
        Child of BiomolComplex. Distributes specified instances across 
        the surface of a sphere. Useful example is the leaflet of a 
        vesicle.
        
    Attributes: 
        Interited
    '''
    def __init__(self, radius, height, number, composition, 
                 pore_radius = 0., seed = None):
        super().__init__(self)
        self.collect_blocks(composition = composition)
        self.generate_positions(radius, height, number)
        if (pore_radius != 0.):
            self.make_equilibration_pores(radius = pore_radius)
        self.gen_ids(seed = seed)
        self.sort_ids()
        self.populate_structure_components()

    def generate_positions(self, radius, height, number):
        '''
        Purpose:   Distribute points across the surface of a sphere.
        Arguments: self) BiomolComplex instance.
                   radius) Float. Radius of the sphere points are to 
                   be distributed on.
                   height) Float. Amount to shift the z-direction.
                   number) Integer. Number of positions to distribute 
                   across the surface of the sphere.
        Requires:  None
        Modifies:  positions
        Returns:   None
        '''
        self.positions = fib_sphere(radius, height, number)

    def make_equilibration_pores(self, axes = ['x', 'y', 'z'], 
                                 radius = 0.0):
        '''
        Purpose:   To make pores for the purpose of equilibrating. The 
                   number of lipids between two vesicle leaflets.
        Arguments: self) BiomolComplex instance.
                   axes) List of Strings. 'x', 'y', and 'z' are 
                   recognized. Determine about which axes pores are 
                   made.
                   radius) Float. Radius of the pores that will be 
                   made.
        Requires:  None
        Modifies:  positions
        Returns:   None
        '''
        if 'x' in axes:
            projection = self.positions[:, 1:]
            proj_mag = np.linalg.norm(projection, axis = 1)
            self.positions = self.positions[proj_mag > radius]
        if 'y' in axes:
            projection = self.positions[:, [0,2]]
            proj_mag = np.linalg.norm(projection, axis = 1)
            self.positions = self.positions[proj_mag > radius]
        if 'z' in axes:
            projection = self.positions[:, :2]
            proj_mag = np.linalg.norm(projection, axis = 1)
            self.positions = self.positions[proj_mag > radius]

    def distribute(self):
        '''
        Purpose:   Orient and place molecules at their asigned 
                   positions. The procedure translates the molecules 
                   to the origin using the tail atom. Aligns the 
                   molecules to their positions using the head atoms. 
                   Distributes the molecules to their positions using 
                   the head atoms.
        Arguments: self) BiomolComplex instance.
                   axes) List of Strings. 'x', 'y', and 'z' are 
                   recognized. Determine about which axes pores are 
                   made.
                   radius) Float. Radius of the pores that will be 
                   made.
        Requires:  positions, entity_atom_count, xyz, head, tail
        Modifies:  xyz, head, tail
        Returns:   None
        '''
        positions = np.repeat(self.positions, self.entity_atom_count, 
                              axis = 0)

        # Translate tail to origin
        self.to_origin()

        #  Rotate head to position
        rot_quat, inv_quat = get_rot_quat(self.head, positions)
        self.rot_quat_mult(rot_quat, inv_quat)

        # Translate head to position
        self.to_positions(positions)
