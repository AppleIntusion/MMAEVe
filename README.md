# MMAEVe - Membranes, Micelles, And Even Vesicles

MMAEVe is a simple python package for creating complex biomolecular systems. It can be used to distribute biomolecules about different surfaces, remove overlap between constructed system elements, write pdb files suitable for use as initial starting structures for AMBER and Gromacs, and write Gromacs topology files. It is a simple but powerful tool that can be used to systematically generate complex structures.

## Installation

Working on getting it up on conda-forge as a package. For now just download the repo. Be sure that your environment has numpy, scipy, and pandas.

## Dirty Guide

Quick and dirty guide. This will be replaced with a much simpler example and this content plus more will be moved to the tutorial.

In the directory where MMAEVE is housed create two directories.
```bash
mkdir structures
mkdir compositions
```
Place any building-blocks you intend to use for making the system in `structures`.

Next we will need to write some composition files. This is by far the most tedious part of using MMAEVe. An example file `compositions/upper_leaf_comp` is found below.

```
POPC     POPS     CHOL     POP2
0.6      0.2      0.12     0.08
2-1-POPC 2-1-POPS 1-1-CHOL 4-1-POP2
7-1-POPC 7-1-POPS 7-1-CHOL 11-1-POP2
```

The top column refers to the names of .pdb files located in the structures directory. The second column should sum to 1. It is the proportion of total molecules that the associated structure should be. The third and fourth lines are of the same format 'atom\_serial\_number-residue\_number-residue\_name. The third and fourth lines, defines which atom should be considered as the "head" and "tail", respectively.

Files in the structures directory
```
$ ls structures
POPC.pdb POPS.pdb POP2.pdb CHOL.pdb c4hre_neu_c4.pdb
```

We will need two more composition files `compositions/lower_leaf_comp` and `compositions/a2t_comp`.

```
POPC     POPS     CHOL     POP2
0.6      0.2      0.12     0.08
7-1-POPC 7-1-POPS 7-1-CHOL 11-1-POP2
2-1-POPC 2-1-POPS 1-1-CHOL 4-1-POP2
```

```
c4hre_neu_c4
1.0
190-90-LEU
1350-956-LYS
```

We can now build out example system.

```python
import MMAEVe as mav

upper_leaf_comp   = mav.read_comp("compositions/upper_leaf_comp")
lower_leaf_comp   = mav.read_comp("compositions/lower_leaf_comp")
a2t_comp          = mav.read_comp("compositions/a2t_comp")
```
The `lower_leaf_comp` inverts the head and tail selection of the `upper_leaf_comp` so that the lipid tails of each bilayer are adjacent to each other.

The three compositions have been imported and the system can now be built. Two instances of `Sphere` are initialized. This will create two leaflets. The first with a radius of 125. Å, height of 0. Å, 3000 lipids, the composition described by `upper_leaf_comp`, and equilibration pores needed for the procedure described by (Reference).

```python
outer_leaf = mav.Sphere(125., 0., 3000, upper_leaf_comp,
                        pore_radius = 20.0)
inner_leaf = mav.Sphere(100., 0., 1875, lower_leaf_comp,
                        pore_radius = 20.0)
```

Lipids are distributed to their positions on the surface of the sphere and the systems are combined to create a vesicle.
```python
outer_leaf.distribute()
inner_leaf.distribute()
vesicle = outer_leaf + inner_leaf
```

We follow a similar procedure to create a flat lipid bilayer but this time we specify the creation of a 350. x 300. Å grid where the upper leaflet is placed 23. Å abover the lower leaflet. Each leaflet has 2008 lipids. The compositions are the same as was used for creating the vesicle.
```python
upper_leaf = mav.Lattice(350., 300., 23., 2008, upper_leaf_comp)
lower_leaf = mav.Lattice(350., 300., 0., 2008, lower_leaf_comp)
upper_leaf.distribute()
lower_leaf.distribute()
bilayer = upper_leaf + lower_leaf
```

The vesicle is translated to to the bilayer then elevated 300. Å above the bilayer.
```python
vesicle + (bilayer.centroid() - vesicle.centroid())
vesicle + np.array([0., 0., 300.])
```

A 200. x 200. Å grid of 9 (3 x 3) A2t is created to place inbetween and, hopefully, tether the vesicle to the bilyaer. We perform some familiar operations.
``` python
a2t = mav.Grid(200., 200., 0., 3, 3, a2t_comp)
a2t.distribute()
a2t + (bilayer.centroid() - a2t.centroid())
a2t + np.array([0., 0., 100.])
```

The three components of the system are combined into a single system.
``` python
vesi_bi_bridge = vesicle + bilayer + a2t
```

The system can be written as either a .cif or .pdb file. The former writes much more quickly. Useful for rapid prototyping of a system. The .cif file will be read by PyMol but not VMD. It has not yet been tested with any other molecular visualization software.
```python
vesi_bi_bridge.write_cif("complexes/vesicle_bi_bridge.cif")
vesi_bi_bridge.write_pdb("complexes/vesicle_bi_bridge.pdb")
```

Will update later. 

## NOTE

This is just a temporary ReadMe.



