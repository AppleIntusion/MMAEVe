'''-----------------------------------------------------*
| Title:   Build Micelle                                |
|                                                       |
| Author:  Gubbin Eel (Satanic Overlord of the Swamp)   |
|                                                       |
| Purpose: Builds a micelle out of lipids.              |
|                                                       |
| NOTE:    Intended improvements are denoted by the #%# |
|          flag. Searching for that flag will identify  |
|          code that needs to be updated or will later  |
|          be improved.                                 |
*-----------------------------------------------------'''

# Import program modules
import manipulate_files as manf
import geom_shapes      as gs

from lipid_structure import lipidStructure

# Import outside modules
import numpy as np
import copy

class Micelle(lipidStructure):
    ''' 
    Description:

        Child class of lipidStructure. The methods defined here
        relate explicitly to building a micelle.

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
    #*bilayer_comp = {"POC" : 0.7, "POS" : 0.1, 
    #*                "PIP" : 0.1, "CHL" : 0.1}
    bilayer_comp = {"POPC" : 0.6,  "POPS" : 0.2, 
                    "CHOL" : 0.12, "POP2" : 0.08}
    # Required information about the proteins present in the system.
    #protein_comp = {"APO" : {"Fraction" : 1.0,
    #                         "Head"     : "29-2-SER",
    #                         "Tail"     : "316-19-ASN"}}
    protein_comp = {"1W7B" : {"Fraction" : 1.0,
                             "Head"     : "7-21-PRO",
                             "Tail"     : "2081-281-LYS"}}
    #protein_comp = {"COV" : {"Fraction" : 1.0,
    #                         "Head"     : "3446-500-THR",
    #                         "Tail"     : "16573-1147-SER"}}
    micelle = Micelle()
    # Add lipid structures and proportions.
    micelle.collect_lipid_blocks(leaf_1 = bilayer_comp)
    # Add protein structures and proportions.
    #*micelle.collect_protein_blocks(proteins = protein_comp)
    # Add lipid positions
    micelle.generate_lipid_points(leaf_1_radius = 150, 
                                  leaf_1_number = 7500, 
                                  structure     = "Micelle")
    # Add protein positions.
    #*micelle.generate_protein_points(leaf_1_radius = 150,
    #*                                leaf_1_number = 10,
    #*                                structure     = "Micelle")
    # Assign IDs to each of the points.
    micelle.gen_ids()
    # Build the micelle using specified parameters.
    micelle.add_lipids()
    # Add proteins to the Micelle system.
    #*micelle.add_proteins()
    # Remove any lipids that overlap with protein.
    #*micelle.remove_overlap()
    # Generate the file lines for the system.
    micelle.struc_to_lines()
    # Save the generate structure to a pdb file.
    micelle.save_lines("test.pdb")
