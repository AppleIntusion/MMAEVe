'''-----------------------------------------------------------*
| Title:   Build Vesicle                                      |
|                                                             |
| Author:  Gubbin Eel (Satanic Overlord of the Swamp)         |
|                                                             |
| Purpose: Builds a vesicle out of lipids.                    |
|                                                             |
| NOTE:    Intended improvements are denoted by the #%# flag. |
|          Searching for that flag will identify code that    |
|          needs to be updated or will later be improved.     |
*-----------------------------------------------------------'''

import manipulate_files as manf
import geom_shapes      as gs
import build_micelle    as bmi

import copy

class Vesicle(bmi.Micelle):
    ''' 
    Inherits the Micelle class. Has two lists as attributes. 
    One is of lipid objects that make up the structure. The 
    other is a list of those lipids as strings in the format 
    of pdb file lines.

    Attributes: struc and lines
    '''

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
            # Points of the sphere to translate Lipid head to.
            sphere_point = self.leaf_1_points[ii]
            # Translate Lipid tail to the origin.
            lipid.trans_axis_to_point([0., 0., 0.], "Tail")
            # Align Lipid head with position vector.
            lipid.align_axis_to_vec(sphere_point, "Head")
            # Translate the lipid to the position on the sphere.
            lipid.trans_axis_to_point(sphere_point, "Head")
            # Add lipid to list
            self.lipids.append(lipid)
            ii += 1

        ii = 0
        for jj in self.leaf_2_ids:
            # Copy Lipid.
            lipid = copy.deepcopy(self.leaf_2_structures[jj])
            # Points of the sphere to translate Lipid head to.
            sphere_point = self.leaf_2_points[ii]
            # Translate Lipid tail to the origin.
            lipid.trans_axis_to_point([0., 0., 0.], "Head")
            # Align Lipid head with position vector.
            lipid.align_axis_to_vec(sphere_point, "Tail")
            # Translate the lipid to the position on the sphere.
            lipid.trans_axis_to_point(sphere_point, "Tail")
            # Add lipid to list
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
            sphere_point = self.protein_points[ii]
            # Translate Protein tail to origin.
            protein.trans_axis_to_point([0., 0., 0.], "Tail")
            # Align Protein head with position vector.
            protein.align_axis_to_vec(sphere_point, "Head")
            # Translate the protein to the position on the sphere.
            protein.trans_axis_to_point(sphere_point, "Head")
            # Add protein to list.
            self.proteins.append(protein)
            ii += 1

if __name__ == "__main__":
    # Information about Lipid bilayer composition and proportions.
    # Initialize Micelle class instance.
    leaf_1_comp = {"POPC" : 0.6,  "POPS" : 0.2, 
                   "CHOL" : 0.12, "POP2" : 0.08}
    leaf_2_comp = {"POPC" : 0.6,  "POPS" : 0.2, 
                   "CHOL" : 0.12, "POP2" : 0.08}
    # Required information about the proteins present in the system.
    #protein_comp = {"APO" : {"Fraction" : 1.0,
    #                         "Head"     : "29-2-SER",
    #                         "Tail"     : "316-19-ASN"}}
    protein_comp = {"1W7B" : {"Fraction" : 1.0,
                             "Head"     : "157-105-PRO",
                             "Tail"     : "559-285-LYS"}}
    #protein_comp = {"COV" : {"Fraction" : 1.0,
    #                         "Head"     : "3446-500-THR",
    #                         "Tail"     : "16573-1147-SER"}}
    vesicle = Vesicle()
    # Add lipid structures and proportions.
    vesicle.collect_lipid_blocks(leaf_1 = leaf_1_comp, 
                                 leaf_2 = leaf_2_comp)
    # Add protein structures and proportions.
    vesicle.collect_protein_blocks(proteins = protein_comp)
    # Add lipid positions
    vesicle.generate_points(leaf_1_radius = 150, 
                            leaf_1_number = 3911, 
                            leaf_2_radius = 125,
                            leaf_2_number = 2716,
                            protein_radius = 240,
                            protein_number = 10,
                            structure     = "Vesicle")
    # Assign IDs to each of the points.
    vesicle.gen_ids()
    # Build the vesicle using specified parameters.
    vesicle.add_lipids()
    # Add proteins to the Micelle system.
    vesicle.add_proteins()
    # Shift the proteins such that they will surround the 
    # vesicle.
    pc = 0
    #for ii in range(len(vesicle.proteins)):
    #    #bilayer.proteins[pc].translate([-12., -12., 175.])
    #    #vesicle.proteins[pc].translate([266.364, 266.300, 266.403])
    #    vesicle.proteins[pc].translate([263.837, 259.993, 262.514])
    #    pc += 1
    # Remove any lipids that overlap with protein.
    #vesicle.remove_overlap()
    # Generate the file lines for the system.
    vesicle.struc_to_lines()
    # Save the generate structure to a pdb file.
    vesicle.save_lines("test.pdb")
