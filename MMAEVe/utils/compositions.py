'''
Compositions
------------

Functions for dealing with compositions. Formatting it as a 
file was much more approachable than having to use a dictionary.
'''

def read_comp(file_name):
    ''' 
    Purpose:   Helper function to read a file and convert it to a 
               composition dictionary. File should follow the 
               following formatting convention:
    
               POPC  POPC  0.6  0.0  2-1-POPC  7-1-POPC 
               POPS  POPS  0.2  0.0  2-1-POPS  7-1-POPS 
               CHOL  CHOL  0.12 0.0  1-1-CHOL  7-1-CHOL 
               POP2  POP2  0.08 0.0  4-1-POP2  11-1-POP2
    
               Where the first column is the prefix of the pdb file 
               name that will be imported. The second column is the 
               name that will be used when writing a gromacs top file. 
               The third column is the fraction of positons that each 
               structure will occupy. The fourth column is the 
               shift relative to the shape normal. The fifth column 
               is the 'serial-residue_number-residue_name' used to 
               specify which atom will serve as the head. The sixth 
               column is the 'serial-residue_number-residue_name' used 
               to specify which atom will serve as the tail.
    Arguments: file_name) String. Name of the composition file.
    Returns:   Dictionary. Properly formatted to be used when 
               initializing any of the supported shapes.
    '''
    # Read file
    with open(file_name, 'r') as f:
        contents = f.readlines()
    # Get fields and info for each individual dict
    comp = {}
    for line in contents:
        struc_data = line.split()
        comp[struc_data[0]] = {"itp"      : struc_data[1], 
                               "Fraction" : float(struc_data[2]), 
                               "Shift"    : float(struc_data[3]), 
                               "Head"     : struc_data[4], 
                               "Tail"     : struc_data[5]}
    return(comp)

