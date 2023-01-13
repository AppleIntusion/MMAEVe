'''-----------------------------------------------------------*
| Title:   Manipulate Files                                   |
|                                                             |
| Author:  Gubbin Eel (Satanic Overlord of the Swamp)         |
|                                                             |
| Purpose: Defines functions used to read, write, and modify  |
|          files. This includes general functions for working |
|          with any file/file-type as well as specific file-  |
|          types such as pdb files.                           |
*-----------------------------------------------------------'''

import molecule_structures as mols
import numpy as np

wspc = "\t\n\x0b\x0c\r "

# Return the contents of a file as a list of strings.
def import_file_contents(file):
    ''' 
    Purpose:   Return the contents of a file as a list of strings.

    Arguments: file) String. The name of a file with its path, if 
               required.
    Returns:   The lines of the file as a list of strings.
    '''
    with open(file, 'r') as f:
        contents = f.readlines()
    return contents

def write_file(file, file_list):
    ''' 
    Purpose:   Write a list of strings to a file.
    Arguments: file) String. The name of the output file with its 
               full path, if required.
               file_list) List of strings. A list containing strings
               representing the lines of each file to be written.
    Returns:   Nothing. Writes a file.
    '''
    with open(file, 'w') as f:
        for line in file_list:
            f.writelines(line)

class PdbFile():
    '''
    Description:
        Used to store all the information found within a PDB file as
        an object that can be manipulated in various ways. Most 
        importantly, the pdb file can be converted to an instance of 
        the Lipid or Protein class.

    Attributes: 
        record_name) PDB record field.
        serial) PDB ATOM record field: Atom serial number.
        atom_name) PDB ATOM record field: Atom name.
        alternate_loc) PDB ATOM record field: Alternate location 
        indicator.
        residue_name) PDB ATOM record field: Residue name.
        chain) PDB ATOM record field: Chain identifier.
        residue_num) PDB ATOM record field: Residue sequence number.
        insertion) PDB ATOM record field: Insertion code.
        x) PDB ATOM record field: x-coordinate.
        y) PDB ATOM record field: y-coordinate.
        z) PDB ATOM record field: z-coordinate.
        occupancy) PDB ATOM record field: Occupancy.
        temp_fac) PDB ATOM record field: Temperature factor.
        element) PDB ATOM record field: Element symbol.
        charge) PDB ATOM record field: Charge.
        line) Original PDB file lines.
    '''

    #%# Update to include a count of file lines and fields. At the
    #   moment I just need it to import ATOM information.
    def __init__(self, file_name = "None", lipid_name = "None", 
                 protein_name = "None"):
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
        ###LIPID MANIPULATION###
        self.record_name   = np.array([])
        self.serial        = np.array([])
        self.atom_name     = np.array([])
        self.alternate_loc = np.array([])
        self.residue_name  = np.array([])
        self.chain         = np.array([])
        self.residue_num   = np.array([])
        self.insertion     = np.array([])
        self.x             = np.array([])
        self.y             = np.array([])
        self.z             = np.array([])
        self.occupancy     = np.array([])
        self.temp_fac      = np.array([])
        self.element       = np.array([])
        self.charge        = np.array([])

        # Assess if a proper number of arguments was provided. Update
        # the decision-making variables accordingly.
        none_count = 0
        if file_name == "None"    : none_count = none_count + 1 
        if lipid_name == "None"   : none_count = none_count + 1 
        if protein_name == "None" : none_count = none_count + 1 
        if none_count != 2:
            print("Incorrect number of arguments specified. Only " + \
                  "specify one of the available optional argument.")
            return
        # Generate file-name if a proper number of arguments was 
        # provided.
        else:
            if file_name != "None":
                pass
            elif lipid_name != "None":
                file_name = "../lipids/" + lipid_name + ".pdb"
            else:
                file_name = "../proteins/" + protein_name + ".pdb" 

        file_lines = import_file_contents(file_name)
        
       
        for line in file_lines:
            if line[0:4] == "ATOM":
                '''
                Atom fields with their index number in the list.
                0:  "ATOM  "        1:  Atom serial num 2:  Atom name 
                3:  Alternate loc   4:  Residue Name    5:  Chain ID 
                6:  Residue seq num 7:  iCode           8:  x-coord 
                9:  y-coord         10: z-coord         11: Occupancy 
                12: Temp factor     13: Elem sym        14: Atomic charge
                '''
                # All line entries for an ATOM entry in a .pdb file.
                line      = [line[0:6],   line[6:11],  line[12:16],
                             line[16],    line[17:21], line[21],
                             line[22:26], line[26],    line[30:38],
                             line[38:46], line[46:54], line[54:60],
                             line[60:66], line[76:78], line[78:80]]

                # Stripping whitespace and formatting numeric entries 
                # to be either integers or floats, whichever is 
                # appropriate.
                self.record_name = np.append(self.record_name, 
                                             line[0].strip(wspc))
                self.serial = np.append(self.serial, 
                                        int(line[1].strip(wspc)))
                self.atom_name = np.append(self.atom_name, 
                                           line[2].strip(wspc))
                self.alternate_loc = np.append(self.alternate_loc, 
                                               line[3].strip(wspc))
                self.residue_name = np.append(self.residue_name, 
                                              line[4].strip(wspc))
                self.chain = np.append(self.chain, 
                                       line[5].strip(wspc))
                self.residue_num = np.append(self.residue_num, 
                                             int(line[6].strip(wspc)))
                self.insertion = np.append(self.insertion, 
                                           line[7].strip(wspc))
                self.x = np.append(self.x, 
                                   float(line[8].strip(wspc)))
                self.y = np.append(self.y, 
                                   float(line[9].strip(wspc)))
                self.z = np.append(self.z, 
                                   float(line[10].strip(wspc)))
                self.occupancy = np.append(self.occupancy, 
                                           float(line[11].strip(wspc)))
                self.temp_fac = np.append(self.temp_fac, 
                                          float(line[12].strip(wspc)))
                self.element = np.append(self.element, 
                                         line[13].strip(wspc))
                try:
                    float(line[14].strip(wspc))
                    self.charge = np.append(self.charge, 
                                            float(line[14].strip(wspc)))
                except ValueError:
                    line[14] = ''
                    self.charge = np.append(self.charge, '')

        self.line = file_lines

    def to_lipid(self):
        ''' 
        Purpose:   Take a PdbFile class instance and convert it to a 
                   Lipid class instance.
        Arguments: self) PdbFile instance.
        Returns:   Lipid instance.
        '''
        x = np.transpose(np.array([self.x]))
        y = np.transpose(np.array([self.y]))
        z = np.transpose(np.array([self.z]))
        xyz = np.concatenate([x, y, z], axis = 1)
        return mols.Lipid(serial = self.serial,       
                          xyz    = xyz, 
                          mass   = np.array([]),      
                          resi   = self.residue_num, 
                          resn   = self.residue_name, 
                          chain  = self.chain, 
                          elem   = self.element,      
                          name   = self.atom_name, 
                          charge = self.charge,       
                          lookup = dict())


    def to_protein(self):
        ''' 
        Purpose:   Take a PdbFile class instance and convert it to a Protein
                   class instance.
        Arguments: self) PdbFile instance.
        Returns:   Protein instance.
        '''
        residues = []
        start = 0
        end   = 0
        current_res_num  = self.residue_num[0]
        previous_res_num = self.residue_num[0]
        for ii in range(len(self.serial)):
            current_res_num = self.residue_num[ii]
            if ((current_res_num != previous_res_num) or 
                  (end + 1 == len(self.serial))):
                if end + 1 == len(self.serial):
                    end += 1
                x = np.transpose(np.array([self.x[start:end]]))
                y = np.transpose(np.array([self.y[start:end]]))
                z = np.transpose(np.array([self.z[start:end]]))
                xyz = np.concatenate([x, y, z], axis = 1)
                res = mols.Residue(serial = self.serial[start:end], 
                                   xyz    = xyz,  
                                   resi   = self.residue_num[start:end], 
                                   resn   = self.residue_name[start:end], 
                                   chain  = self.chain[start:end], 
                                   elem   = self.element[start:end], 
                                   name   = self.atom_name[start:end], 
                                   charge = self.charge[start:end])
                residues.append(res)
                start = end

            end += 1

            previous_res_num = current_res_num

        return mols.Protein(residues = np.array(residues))

if __name__ == "__main__":
    #pro = PdbFile(protein_name = "apoa1_peptide")
    pro = PdbFile(protein_name = "APO")
    #pro = pro.to_protein()
    pro = pro.to_protein()
    #print(pro.residues[0].resn)
    write_file("test.pdb", pro.to_pdb_file())
    #for res in pro.residues:
    #    for atom in res:
    #        print(res[atom].serial, ' ', res[atom].resi)
