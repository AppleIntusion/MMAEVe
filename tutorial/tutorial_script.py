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

'''
| Nanodisc
'''
nanodisc_comp = mav.read_comp("compositions/nanodisc_comp")

upper_leaf = mav.Disc(95., 23., 464, upper_leaf_comp)
lower_leaf = mav.Disc(95.,  0., 464, lower_leaf_comp)

upper_leaf.distribute()
lower_leaf.distribute()

disc = upper_leaf + lower_leaf

msp2n2 = mav.Lattice(115., 115., 0., 1, nanodisc_comp)
msp2n2 + (disc.centroid() - msp2n2.centroid())

disc.remove_overlap(msp2n2, radius = 3.0)
nanodisc = disc + msp2n2 

nanodisc.write_cif("complexes/nanodisc.cif")
nanodisc.write_pdb("complexes/nanodisc.pdb")

'''
| Nanotube
'''

outer_leaf = mav.Cylinder(75., 300., 0., 2207, upper_leaf_comp)
inner_leaf = mav.Cylinder(52., 300., 0., 1606, lower_leaf_comp)

outer_leaf.distribute()
inner_leaf.distribute()

nanotube = outer_leaf + inner_leaf

nanotube.write_cif("complexes/nanotube.cif")
nanotube.write_pdb("complexes/nanotube.pdb")

'''
| Array of Nanotubes
'''

nanotube0 = copy.deepcopy(nanotube)
nanotube1 = copy.deepcopy(nanotube)
nanotube2 = copy.deepcopy(nanotube)

nanotube0 + np.array([165., 0., 0.])
nanotube1 + np.array([0., 165., 0.])
nanotube2 + np.array([165., 165., 0.])

nanotube_array = nanotube + nanotube0 + nanotube1 + nanotube2

nanotube_array.write_cif("complexes/nanotube_array.cif")
nanotube_array.write_pdb("complexes/nanotube_array.pdb")

'''
| Vesicle
'''

outer_leaf = mav.Sphere(125., 0., 3000, upper_leaf_comp, 
                    pore_radius = 20.0)
inner_leaf = mav.Sphere(100., 0., 1875, lower_leaf_comp,
                    pore_radius = 20.0)

outer_leaf.distribute()
inner_leaf.distribute()

vesicle = outer_leaf + inner_leaf

vesicle.write_cif("complexes/vesicle.cif")
vesicle.write_pdb("complexes/vesicle.pdb")

'''
| Periphreal Membrane-Binding Proteins around a Vesicle
'''

a2_comp = mav.read_comp("compositions/a2_comp")

a2 = mav.Sphere(165., 0., 10, a2_comp)
a2.distribute()

vesi_a2 = vesicle + a2

vesi_a2.write_cif("complexes/vesi_a2.cif")
vesi_a2.write_pdb("complexes/vesi_a2.pdb")

'''
| Membrane-Vesicle Junction
'''

a2t_comp = mav.read_comp("compositions/a2t_comp")

vesicle0 = copy.deepcopy(vesicle)

upper_leaf = mav.Lattice(350., 350., 23., 2008, upper_leaf_comp)
lower_leaf = mav.Lattice(350., 350., 0., 2008, lower_leaf_comp)
upper_leaf.distribute()
lower_leaf.distribute()
bilayer = upper_leaf + lower_leaf

vesicle0 + (bilayer.centroid() - vesicle0.centroid())
vesicle0 + np.array([0., 0., 300.])

a2t = mav.Grid(240., 240., 0., 3, 3, a2t_comp)
a2t.distribute() 
a2t + (bilayer.centroid() - a2t.centroid())
a2t + np.array([0., 0., 100.])

vesi_bi_a2t = vesicle0 + bilayer + a2t

vesi_bi_a2t.write_cif("complexes/vesicle_bi_a2t.cif")
vesi_bi_a2t.write_pdb("complexes/vesicle_bi_a2t.pdb")

'''
| Covid Virion
'''

spike_comp = mav.read_comp("compositions/spike_comp")

outer_leaf = mav.Sphere(500., 0., 51475, upper_leaf_comp)
inner_leaf = mav.Sphere(475., 0., 46456, lower_leaf_comp)
outer_leaf.distribute()
inner_leaf.distribute()
vesicle = outer_leaf + inner_leaf

spike = mav.Sphere(600., 0., 75, spike_comp)
spike.distribute()

vesicle.remove_overlap(spike, 4.)

covid = spike + vesicle
covid.write_cif("complexes/covid_viron.cif")
covid.write_pdb("complexes/covid_viron.pdb")

'''
| GROMACS Topology
'''

bilayer.write_gromacs_top("complexes/bilayer.top")
