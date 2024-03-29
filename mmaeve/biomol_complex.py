'''
Biomolecular Complex

Code related to building biomolecular systems.
'''

import numpy as np
from .pdb_file import *
from .utils import *
import decimal as dec
import scipy

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
                 ratios = None, itp_names = None, shifts = None, 
                 ids = None, serial = None, name = None, resi = None, 
                 resn = None, chain = None, xyz = None, elem = None, 
                 mass = None, charge = None, radius = None, 
                 head = None, tail = None, entity_count = 0, 
                 entity_atom_count = None, residue_count = None, 
                 residue_atom_count = None):
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
        if itp_names is None:
            self.itp_names = np.array([]) 
        else:
            self.itp_names = itp_names 
        if shifts is None:
            self.shifts = np.array([]) 
        else:
            self.shifts = shifts
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
        itp_names          = np.concatenate([self.itp_names, 
                                             struc2.itp_names])
        shifts             = np.concatenate([self.shifts, 
                                             struc2.shifts])
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
                structures = structures, ratios = ratios, 
                itp_names = itp_names, ids = ids, serial = serial, 
                name = name, resi = resi, resn = resn, chain = chain, 
                xyz = xyz, elem = elem, mass = mass, charge = charge, 
                radius = radius, head = head, tail = tail, 
                entity_count = entity_count, 
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
        self.itp_names = []
        self.shifts = []
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
            self.itp_names.append(composition[name]["itp"])
            self.shifts.append(composition[name]["Shift"])
        self.structures = np.array(self.structures)
        self.ratios = np.array(self.ratios)
        self.itp_names = np.array(self.itp_names)
        self.shifts = np.array(self.shifts)
        # Check for instances of 0% composition and remove those 
        # structures.
        self.itp_names = self.itp_names[self.ratios != 0.]
        self.structures = self.structures[self.ratios != 0.]
        self.ratios = self.ratios[self.ratios != 0.]

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

        total_lipids = len(self.positions)
        id_pool = np.array([])

        dec.getcontext().prec = 10
        proportions = []
        numbers = []
        for ii, ratio in enumerate(self.ratios):
            proportion = ratio / np.sum(self.ratios)
            proportion = np.round(proportion, 5)
            number = int(np.round(proportion * total_lipids, 0))
            id_pool = np.concatenate([id_pool, np.repeat(ii, number)])
        id_pool = id_pool.astype(int)

        # It pains me to admit it, but the rounding error spanked me.
        # This is a temporary fix, if there are a mismatch of positions
        # to ids then the number is corrected. Shouldn't ever be more
        # than one but, the solution will handle cases where the total
        # number is off by more than one. DEBUG
        if len(id_pool) < len(self.positions):
            to_remove = len(self.positions) - len(id_pool)
            for ii in range(to_remove):
                self.positions = np.delete(self.positions, -1, 
                                           axis = 0)
        elif len(id_pool) > len(self.positions):
            to_remove = len(id_pool) - len(self.positions)
            for ii in range(to_remove):
                id_pool = np.delete(id_pool, -1, axis = 0)
        
        if len(id_pool) != len(self.positions):
            print("YOUR CODE IS BROKEN. CHECK `gen_ids`", len(id_pool), len(self.positions)) # DEBUG

        state.shuffle(id_pool)
        self.ids = id_pool

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
            chunks = list(range(0, chunksize, nparts))
            chunks.append(nparts + 1)   # so slicing is from 2ndtolast:size

            for startind, stopind in zip(chunks[:-1], chunks[1:]):
                fout.writelines(
                    [ii[9].format(ii[0], ii[1], ii[2], ii[3], 
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
        itp_names = self.itp_names
        _, counts = np.unique(self.ids, return_counts = True)
        counts = counts.astype(str)
        top_lines = np.char.add(itp_names, '\t\t')
        top_lines = np.char.add(top_lines, counts)
        top_header = ["; You'll need to add all of the itp-file " + \
                      "references", " ", "[ molecules ]", 
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
        self.apply_shift()
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

    def apply_shift(self):
        '''
        Purpose:   Apply the previously specified shift to the system.
        Arguments: self) BiomolComplex instance.
        Requires:  None
        Modifies:  positions
        Returns:   None
        '''
        self.positions[:, 2] = self.positions[:, 2] + \
                               self.shifts[self.ids]
        

    def distribute(self):
        '''
        Purpose:   Orient and place molecules at their asigned 
                   positions. The procedure translates the molecules 
                   to the origin using the tail atom. Aligns the 
                   molecules to the z-axis. Distributes the molecules 
                   to their positions using the head atoms.
                   itp_names = self.itp_names
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
        self.apply_shift()
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
        self.apply_shift()
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
        self.apply_shift()
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
        self.apply_shift()
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

    def apply_shift(self):
        '''
        Purpose:   Apply the previously specified shift to the system.
        Arguments: self) BiomolComplex instance.
        Requires:  None
        Modifies:  positions
        Returns:   None
        '''
        shift = np.copy(self.positions)
        shift[:, 2] = 0.
        mag   = np.linalg.norm(shift, axis = 1).reshape((len(self.positions), 1))
        shift = shift / mag
        shift = shift * self.shifts[self.ids].reshape((len(self.positions), 1))       
        self.positions = self.positions + shift

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
        self.apply_shift()
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

    def apply_shift(self):
        '''
        Purpose:   Apply the previously specified shift to the system.
        Arguments: self) BiomolComplex instance.
        Requires:  None
        Modifies:  positions
        Returns:   None
        '''
        shift = np.linalg.norm(self.positions, axis = 1).reshape((len(self.positions), 1))
        shift = self.positions / shift
        shift = shift * self.shifts[self.ids].reshape((len(self.positions), 1))       
        self.positions = self.positions + shift

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
