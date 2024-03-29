'''
Molecules

Code related to the representation and manipulation of molecules.
'''

import numpy as np

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
