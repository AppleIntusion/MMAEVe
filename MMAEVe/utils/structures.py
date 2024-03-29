'''
Structures
----------

Helper functions for dealing with structures that fall outside of 
the scope of the main organization of MMAEVe. This is primarily
to deal with Martini3 lipids. PDBs are no longer available via
their website but are incorporated into the Insane script.
These helper functions allow for one to convert one of the Insane
templates into a pdb file for use with MMAEVe.
'''

def insane_to_pdb(template):
    ''' 
    Purpose:   Helper function to read a single dictionary 
               approximately following the template style of insane 
               lipids. "Approximately" refers to the fact that the 
               style is not exactly that used by Insane but it is 
               quite easy to make the dictionary in a text editor if 
               the insane template is available. Only handles one 
               lipid at a time. An example follows.
    
               The original template should look something like this:

               moltype = "sterol"
               lipidsx[moltype] = (0, 1, 0, 0, 1, 0, 0.5, 0.5, 0, 0)
               lipidsy[moltype] = (0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
               lipidsz[moltype] = (5.3, 4.5, 3.9, 3.3, 3, 
                                   2.6, 4.5, 2.6, 1.4, 0)
               template = " ROH  R1  R2  R3  R4   - R5  R6  C1  C2 "
               lipidsa.update({"CHOL": (moltype, template),

               You will need to convert it to the following form:

               template = {'x' : (0, 1, 0, 0, 1, 0, 0.5, 0.5, 0, 0),
                           'y' : (0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                           'z' : (5.3, 4.5, 3.9, 3.3, 3, 2.6, 4.5, 2.6, 1.4, 0),
                           "lip_temp" : "ROH R1 R2 R3 R4 - R5 R6 C1 C2"
                           "lip_name" : "CHOL"
               }

               where the "lip_name" entry will be used as the residue
               name in the .pdb file.
    Arguments: template) Dictionary. Should be formatted as discussed
               in `Purpose`.
    Returns:   Dictionary. Contains all info needed to write a pdb of
               the template.
    '''
    atom_names = template["lip_temp"].split()
    x_coords   = [template['x'][ii] for ii in range(len(atom_names)) \
                  if atom_names[ii] != '-']
    y_coords   = [template['y'][ii] for ii in range(len(atom_names)) \
                  if atom_names[ii] != '-']
    z_coords   = [template['z'][ii] for ii in range(len(atom_names)) \
                  if atom_names[ii] != '-']
    atom_names = [atom_names[ii] for ii in range(len(atom_names)) \
                  if atom_names[ii] != '-']
    resn = np.repeat(template["lip_name"], len(x_coords))
    pdb_data = {"resn"  : resn,
                "names" : atom_names,
                'x'     : x_coords,
                'y'     : y_coords,
                'z'     : z_coords}
    return(pdb_data)

def write_pdb(pdb_data, file_name, buff = 8192, chunksize = 100000):
    '''
    Purpose:   Take the output of `insane_to_pdb` and write a .pdb
               file.
    Arguments: pdb_data) Dictionary. Output of `insane_to_pdb`.
               file_name) String. Name to write the pdb file to.
               should include the '.pdb' extension in the name.
               buff) Integer. Controls the buffersize of the open
               function. 
               chunksize) Integer. Controls the number of "chunks" 
               to write the file in. 
    Returns:   Nothing
    '''
    atom_count = len(pdb_data['x'])            # Total number of atoms
    residue_atom_count = atom_count
    res_count = 1 # Total # of residues
    write_string = "{:4s}  {:5d} {:4s} {:4s}{:1s}{:4d}    " + \
                   "{:8.3f}{:8.3f}{:8.3f}\n"

    # PDB fields to write: atom serial, atom name, residue name, 
    # chain id, residue number, x, y, z.
    fields = [
        np.repeat("ATOM", atom_count),
        np.arange(1, atom_count + 1, 1), # len correctable
        pdb_data["names"],
        pdb_data["resn"],
        np.repeat('', len(pdb_data["resn"])),
        np.repeat(np.arange(1, res_count + 1, 1), 
                  residue_atom_count), # len correctable
        np.round(pdb_data['x'], decimals = 3),
        np.round(pdb_data['y'], decimals = 3),
        np.round(pdb_data['z'], decimals = 3),
        np.repeat(write_string, atom_count)]
    # Reduce atom serial and residue number to proper size
    fields[1] = np.mod(fields[1], 100000)
    fields[5] = np.mod(fields[5], 10000)

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
