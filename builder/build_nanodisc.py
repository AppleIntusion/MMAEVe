'''-----------------------------------------------------------*
| Title:     Build Nanodisc                                   |
|                                                             |
| Author(s): Dumpster Monkey (Gaurdian of the Trailer Park)   |
|            Gubbin Eel (Satanic Overlord of the Swamp)       |
|                                                             |
| Purpose:   Builds a nanodisc out of lipids an small         |
|            segments of scaffold protein.                    |
|                                                             |
| NOTE:      Intended improvements are denoted by the #%#     |
|            flag. Searching for that flag will identify code |
|            that needs to be updated or will later be        |
|            improved.                                        |
*-----------------------------------------------------------'''

# Import program modules
import manipulate_files as manf
import geom_shapes      as gs
import lipid_structure as ls

# Import outside modules
import numpy as np
import copy

class Nanodisc(ls.lipidStructure):
    ''' 
    Description:

        Child class of lipidStructure. The methods defined here
        relate explicitly to building a nandisc.

    Attributes: 
 
        See lipidStructure class attributes.
    '''
    
    def __init__(self):
        super().__init__()

    def add_lipids(self):
        '''
        Purpose:
        Arguments:
        Returns:
        '''
        ii = 0
        for jj in self.leaf_1_ids:
            # Copy Lipid.
            lipid = copy.deepcopy(self.leaf_1_structures[jj])
            # Points of the sphere to translate Lipid head to
            plane_point = self.leaf_1_points[ii]
            # Translate Lipid tail to the origin
            lipid.trans_axis_to_point([0., 0., 0.], "Tail")
            # Align vector to known vector
            lipid.align_axis_to_vec([1., 1., 0], "Head")
            # Rotate Lipid such that it is aligned to [0., 0., 1.]
            lipid.rotate(axis = 'z', theta = 3.14 * 0.25)
            lipid.rotate(axis = 'x', theta = 3.14 * 0.5)
            # Translate to position on the grid
            lipid.trans_axis_to_point(plane_point, "Head")
            self.lipids.append(lipid)
            ii += 1

        ii = 0
        for jj in self.leaf_2_ids:
            # Copy Lipid.
            lipid = copy.deepcopy(self.leaf_2_structures[jj])
            # Points of the sphere to translate Lipid head to
            plane_point = self.leaf_2_points[ii]
            # Translate Lipid tail to the origin
            lipid.trans_axis_to_point([0., 0., 0.], "Tail")
            # Align vector to known vector
            lipid.align_axis_to_vec([1., 1., 0], "Head")
            # Rotate Lipid such that it is aligned to [0., 0., 1.]
            lipid.rotate(axis = 'z', theta = 3.14 * 0.25)
            lipid.rotate(axis = 'x', theta = 3.14 * 1.5)
            # Translate to position on the grid
            lipid.trans_axis_to_point(plane_point, "Head")
            self.lipids.append(lipid)
            ii += 1

    def add_proteins(self):
        '''
        Purpose:
        Arguments:
        Usage:
        '''
        ii = 0
        for jj in self.protein_ids:
            # Copy Protein.
            protein = copy.deepcopy(self.protein_structures[jj])
            # Points of the sphere to translate Protein head to.
            plane_point = self.protein_points[ii]
            # Translate Protein tail to origin.
            protein.trans_axis_to_point([0., 0., 0.], "Tail")
            # Align vector to known vector
            protein.align_axis_to_vec([1., 1., 0], "Head")
            protein.align_axis_to_vec([1., 1., 0], "Head")
            # Rotate Lipid such that it is aligned to [0., 0., 1.]
            protein.rotate(axis = 'z', theta = 3.14 * 0.25)
            protein.rotate(axis = 'x', theta = 3.14 * 0.5)
            # Translate to position on the grid
            protein.trans_axis_to_point(plane_point, "Head")
            self.proteins.append(protein)
            ii += 1

if __name__ == "__main__":
    
    nanodisc = Nanodisc()

    '''
    Part I: Build the Nanodisc with the specified concentration.
    '''
    # Information about Lipid bilayer composition and proportions.
    leaf_1_comp = {"POPC" : 0.4,  "POPS" : 0.2, 
                   "CHOL" : 0.12, "POP2" : 0.28}
    leaf_2_comp = {"POPC" : 0.4,  "POPS" : 0.2, 
                   "CHOL" : 0.12, "POP2" : 0.28}
    protein_comp = {"APO" : {"Fraction" : 1.0,
                             "Head"     : "29-2-SER",
                             "Tail"     : "313-19-ASN"}}

    nanodisc = Nanodisc()

    nanodisc.collect_lipid_blocks(leaf_1 = leaf_1_comp, 
                                  leaf_2 = leaf_2_comp)
    nanodisc.collect_protein_blocks(proteins = protein_comp)

    nanodisc.generate_lipid_points(leaf_1_radius = 150., 
                                   leaf_1_number = 1100, 
                                   leaf_2_radius = 150., 
                                   leaf_2_number = 1100,
                                   height        = 43.,  
                                   structure     = "Nanodisc")
    nanodisc.generate_protein_points(protein_radius = 158., 
                                     protein_number = 60, 
                                     pro_height     = 38.,
                                     structure      = "Circle")

    nanodisc.gen_ids()

    nanodisc.add_lipids()
    nanodisc.add_proteins()

    '''
    Part II: Save Structure
    '''
    # Convert points to pdb file lines.
    nanodisc.struc_to_lines()
    # Save the generate structure to a pdb file.
    nanodisc.save_lines("test.pdb")
