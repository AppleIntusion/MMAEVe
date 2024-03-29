'''
PDB File Import

Code related to the import, representation, and conversion of PDB 
files.
'''

import pandas as pd
import numpy as np
from .molecules import *

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
