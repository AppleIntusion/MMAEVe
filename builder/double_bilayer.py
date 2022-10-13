from build_bilayer import *

'''
Part I: Build the Bilayer with Porins Inserted
'''
# Information about Lipid bilayer composition and proportions.
leaf_1_comp = {"POPC" : 0.6,  "POPS" : 0.2, 
               "CHOL" : 0.12, "POP2" : 0.08}
leaf_2_comp = {"POPC" : 0.6,  "POPS" : 0.2, 
               "CHOL" : 0.12, "POP2" : 0.08}
protein_comp = {"2POR" : {"Fraction" : 1.0,
                          "Head"     : "323-158-PHE",
                          "Tail"     : "521-256-THR"}}
bilayer = Bilayer()
bilayer.collect_lipid_blocks(leaf_1 = leaf_1_comp, 
                             leaf_2 = leaf_2_comp)
bilayer.collect_protein_blocks(proteins = protein_comp)
bilayer.generate_points(width  = 350.,
                        length = 350.,
                        height = 43.,
                        leaf_1_number = 1915,
                        leaf_2_number = 1915,
                        pro_width = 250,
                        pro_length = 250,
                        pro_height = 43.,
                        protein_number = 4,
                        structure = "Bilayer")
bilayer.gen_ids()
bilayer.add_lipids()
bilayer.add_proteins()
pc = 0
for ii in range(len(bilayer.proteins)):
    #bilayer.proteins[pc].translate([-12., -12., 175.])
    bilayer.proteins[pc].translate([0., -16., 175.])
    pc += 1
lc = 0
for ii in range(len(bilayer.lipids)):
    bilayer.lipids[lc].translate([0., 0., 175.])
    lc += 1

'''
Part III: Add Another Bilayer Layer with Porins
'''
bilayer.leaf_1_points      = np.array([]) 
bilayer.leaf_2_points      = np.array([]) 
bilayer.protein_points     = np.array([]) 
bilayer.leaf_1_structures  = []
bilayer.leaf_2_structures  = [] 
bilayer.protein_structures = [] 
bilayer.leaf_1_ratios      = np.array([]) 
bilayer.leaf_2_ratios      = np.array([]) 
bilayer.protein_ratios     = np.array([]) 
bilayer.leaf_1_ids         = []
bilayer.leaf_2_ids         = []
bilayer.protein_ids        = []
bilayer.lines              = []
bilayer.collect_lipid_blocks(leaf_1 = leaf_1_comp, 
                             leaf_2 = leaf_2_comp)
bilayer.collect_protein_blocks(proteins = protein_comp)
bilayer.generate_points(width  = 350.,
                        length = 350.,
                        height = 43.,
                        leaf_1_number = 1915,
                        leaf_2_number = 1915,
                        pro_width = 250,
                        pro_length = 250,
                        pro_height = 43.,
                        protein_number = 4,
                        structure = "Bilayer")
bilayer.gen_ids()
bilayer.add_lipids()
bilayer.add_proteins()
for ii in range(len(bilayer.proteins) - pc):
    #bilayer.proteins[pc].translate([-12., -12., 0.])
    bilayer.proteins[pc].translate([0., -16., 0.])
    pc += 1

bilayer.remove_overlap()

'''
Part II: Add 4HRE Tetramers
'''
protein_comp = {"4HRE" : {"Fraction" : 1.0,
                          "Head"     : "42-19-PRO",
                          "Tail"     : "1936-2560-GLN"}}


bilayer.protein_structures = []
bilayer.collect_protein_blocks(proteins = protein_comp)
bilayer.generate_points(width  = 350.,
                        length = 350.,
                        height = 83.,
                        leaf_1_number = 0,
                        leaf_2_number = 0,
                        pro_width = 100,
                        pro_length = 100,
                        pro_height = 80.,
                        protein_number = 4,
                        structure = "Bilayer")
bilayer.add_proteins()
for ii in range(len(bilayer.proteins) - pc):
    bilayer.proteins[pc].translate([0., 0., 25.])
    pc += 1


'''
Part IV: Save Structure
'''
# Convert points to pdb file lines.
bilayer.struc_to_lines()
# Save the generate structure to a pdb file.
bilayer.save_lines("test.pdb")
