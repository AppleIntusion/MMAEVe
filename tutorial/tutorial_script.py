import MMAEVe as mav
import numpy as np
import copy 

'''
| Simple Bilayer
'''

upper_leaf_comp   = mav.read_comp("compositions/upper_leaf_comp")
lower_leaf_comp   = mav.read_comp("compositions/lower_leaf_comp")

upper_leaf = mav.Lattice(200., 200., 23., 655, upper_leaf_comp)
lower_leaf = mav.Lattice(200., 200., 0., 655, lower_leaf_comp)

upper_leaf.distribute()
lower_leaf.distribute()

bilayer = upper_leaf + lower_leaf

bilayer.write_cif("complexes/bilayer.cif")
bilayer.write_pdb("complexes/bilayer.pdb")

'''
| Single Embedded Porin
'''
porin_comp = mav.read_comp("compositions/porin_comp")

bilayer0 = copy.deepcopy(bilayer)

porin = mav.Lattice(200., 200., 0., 1, porin_comp)
porin.distribute() 

porin + (bilayer0.centroid() - porin.centroid())

bilayer0.remove_overlap(porin, 20.0, atom_based = False)
bilayer0.remove_overlap(porin, 3.0, atom_based = True)

porin_bilayer = bilayer0 + porin

porin_bilayer.write_cif("complexes/porin_bilayer.cif")
porin_bilayer.write_pdb("complexes/porin_bilayer.pdb")

'''
| Multipe Embedded Porin
'''

bilayer1 = copy.deepcopy(bilayer)

porin_grid = mav.Grid(150., 150., 0., 2, 2, porin_comp)
porin_grid.distribute() 

porin_grid + (bilayer1.centroid() - porin_grid.centroid())

bilayer1.remove_overlap(porin_grid, 20.0, atom_based = False)
bilayer1.remove_overlap(porin_grid, 3.0, atom_based = True)

porin_grid_bilayer = bilayer1 + porin_grid

porin_grid_bilayer.write_cif("complexes/porin_grid_bilayer.cif")
porin_grid_bilayer.write_pdb("complexes/porin_grid_bilayer.pdb")

### Multiple Proteins








## Lipid Nanodiscs

#nanodisc_comp     = mav.read_comp("compositions/nanodisc_comp")
#spike_comp        = mav.read_comp("compositions/spike_comp")
#a2t_comp          = mav.read_comp("compositions/a2t_comp")
#spike_coarse_comp = mav.read_comp("compositions/cspike_comp")

## Lipid Nanotubes

## An Array of Lipid Nanotubes

## Vesicles

## Periphreal Membrane-Binding Proteins around a Vesicle

## Membrane-Vesicle Junction

## Covid Virion

## GROMACS Topology Files
#bilayer.write_gromacs_top("complexes/bilayer.top")

## AMBER-Safe PDB Files

## Reproducible Systems
