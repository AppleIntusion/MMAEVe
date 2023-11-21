import MMAEVe as mav
import numpy as np
import copy 
import time

def generate_timed(func):
    def timed_func(*args, **kwargs):
        bench_times = []
        for ii in range(5):
            start_time = time.time()
            func(*args, **kwargs)
            bench_times.append(time.time() - start_time)
        print("Mean: %.3fs" % np.mean(bench_times))
        print("Standard Deviation: %.3fs" % np.std(bench_times))
        print('')
    return timed_func

'''
| Simple Bilayer
'''
@generate_timed
def bilayer_benchmark():
    upper_leaf_comp   = mav.read_comp("compositions/upper_leaf_comp")
    lower_leaf_comp   = mav.read_comp("compositions/lower_leaf_comp")
    
    upper_leaf = mav.Lattice(200., 200., 23., 655, upper_leaf_comp, seed = 0)
    lower_leaf = mav.Lattice(200., 200., 0., 655, lower_leaf_comp, seed = 0)
    
    upper_leaf.distribute()
    lower_leaf.distribute()
    
    bilayer = upper_leaf + lower_leaf
    
    bilayer.write_pdb("complexes/bilayer.pdb")

'''
| Multipe Embedded Porin
'''
@generate_timed
def porin_gird_benchmark():
    upper_leaf_comp   = mav.read_comp("compositions/upper_leaf_comp")
    lower_leaf_comp   = mav.read_comp("compositions/lower_leaf_comp")
    
    upper_leaf = mav.Lattice(200., 200., 23., 655, upper_leaf_comp, seed = 0)
    lower_leaf = mav.Lattice(200., 200., 0., 655, lower_leaf_comp, seed = 0)
    
    upper_leaf.distribute()
    lower_leaf.distribute()
    
    bilayer = upper_leaf + lower_leaf

    porin_comp = mav.read_comp("compositions/porin_comp")
    
    porin_grid = mav.Grid(150., 150., 0., 2, 2, porin_comp, seed = 0)
    porin_grid.distribute() 
    
    porin_grid + (bilayer.centroid() - porin_grid.centroid())
    
    bilayer.remove_overlap(porin_grid, 20.0, atom_based = False)
    bilayer.remove_overlap(porin_grid, 3.0, atom_based = True)
    
    porin_grid_bilayer = bilayer + porin_grid
    
    porin_grid_bilayer.write_pdb("complexes/porin_grid_bilayer.pdb")

'''
| Nanodisc
'''
@generate_timed
def nanodisc_benchmark():
    upper_leaf_comp   = mav.read_comp("compositions/upper_leaf_comp")
    lower_leaf_comp   = mav.read_comp("compositions/lower_leaf_comp")

    nanodisc_comp = mav.read_comp("compositions/nanodisc_comp")
    
    upper_leaf = mav.Disc(95., 23., 464, upper_leaf_comp, seed = 0)
    lower_leaf = mav.Disc(95.,  0., 464, lower_leaf_comp, seed = 0)
    
    upper_leaf.distribute()
    lower_leaf.distribute()
    
    disc = upper_leaf + lower_leaf
    
    msp2n2_1 = mav.Lattice(115., 115., 0., 1, nanodisc_comp, seed = 0)
    msp2n2_1 + (disc.centroid() - msp2n2_1.centroid())

    disc.remove_overlap(msp2n2_1, radius = 3.0)
    nanodisc = disc + msp2n2_1
    
    nanodisc.write_pdb("complexes/nanodisc.pdb")

'''
| Nanodisc Spike
'''
@generate_timed
def nanodisc_spike_benchmark():
    upper_leaf_comp   = mav.read_comp("compositions/upper_leaf_comp")
    lower_leaf_comp   = mav.read_comp("compositions/lower_leaf_comp")

    nanodisc_comp = mav.read_comp("compositions/nanodisc_comp")

    spike_comp = mav.read_comp("compositions/spike_comp")
    
    upper_leaf = mav.Disc(95., 23., 464, upper_leaf_comp, seed = 0)
    lower_leaf = mav.Disc(95.,  0., 464, lower_leaf_comp, seed = 0)
    
    upper_leaf.distribute()
    lower_leaf.distribute()
    
    disc = upper_leaf + lower_leaf
    
    msp2n2_1 = mav.Lattice(115., 115., 0., 1, nanodisc_comp, seed = 0)
    msp2n2_1 + (disc.centroid() - msp2n2_1.centroid())

    spike = mav.Lattice(95., 95., 0., 1, spike_comp, seed = 0)
    spike + (disc.centroid() - spike.centroid())
    spike + np.array([0., 0., 70.])

    disc.remove_overlap(msp2n2_1, radius = 3.0)
    disc.remove_overlap(spike, radius = 3.0)
    nanodisc = disc + msp2n2_1 + spike
    
    nanodisc.write_pdb("complexes/nanodisc_spike.pdb")

'''
| Nanotube
'''
@generate_timed
def nanotube_benchmark():
    upper_leaf_comp   = mav.read_comp("compositions/upper_leaf_comp")
    lower_leaf_comp   = mav.read_comp("compositions/lower_leaf_comp")
    
    outer_leaf = mav.Cylinder(50., 250., 0., 1839, upper_leaf_comp, seed = 0)
    inner_leaf = mav.Cylinder(27., 250., 0., 993, lower_leaf_comp, seed = 0)
    
    outer_leaf.distribute()
    inner_leaf.distribute()
    
    nanotube = outer_leaf + inner_leaf
    
    nanotube.write_pdb("complexes/nanotube.pdb")

'''
| Array of Nanotubes
'''

@generate_timed
def nanotube_array_benchmark():
    upper_leaf_comp   = mav.read_comp("compositions/upper_leaf_comp")
    lower_leaf_comp   = mav.read_comp("compositions/lower_leaf_comp")
    
    outer_leaf = mav.Cylinder(50., 250., 0., 1839, upper_leaf_comp, seed = 0)
    inner_leaf = mav.Cylinder(27., 250., 0., 993, lower_leaf_comp, seed = 0)
    
    outer_leaf.distribute()
    inner_leaf.distribute()
    
    nanotube = outer_leaf + inner_leaf

    nanotube0 = copy.deepcopy(nanotube)
    nanotube1 = copy.deepcopy(nanotube)
    nanotube2 = copy.deepcopy(nanotube)
    
    nanotube0 + np.array([120., 0., 0.])
    nanotube1 + np.array([0., 120., 0.])
    nanotube2 + np.array([120., 120., 0.])
    
    nanotube_array = nanotube + nanotube0 + nanotube1 + nanotube2
    
    nanotube_array.write_pdb("complexes/nanotube_array.pdb")

'''
| Vesicle
'''

@generate_timed
def vesicle_benchmark():
    upper_leaf_comp   = mav.read_comp("compositions/upper_leaf_comp")
    lower_leaf_comp   = mav.read_comp("compositions/lower_leaf_comp")
    
    outer_leaf = mav.Sphere(125., 0., 3626, upper_leaf_comp, 
                        pore_radius = 20.0)
    inner_leaf = mav.Sphere(100., 0., 1857, lower_leaf_comp,
                        pore_radius = 20.0)
    
    outer_leaf.distribute()
    inner_leaf.distribute()
    
    vesicle = outer_leaf + inner_leaf
    
    vesicle.write_pdb("complexes/vesicle.pdb")

'''
| Membrane-Vesicle Junction
'''
@generate_timed
def mem_vesi_benchmark():
    upper_leaf_comp   = mav.read_comp("compositions/upper_leaf_comp")
    lower_leaf_comp   = mav.read_comp("compositions/lower_leaf_comp")
    a2t_comp = mav.read_comp("compositions/a2t_comp")
    
    outer_leaf = mav.Sphere(125., 0., 3626, upper_leaf_comp, 
                        pore_radius = 20.0)
    inner_leaf = mav.Sphere(100., 0., 1857, lower_leaf_comp,
                        pore_radius = 20.0)
    
    outer_leaf.distribute()
    inner_leaf.distribute()
    
    vesicle = outer_leaf + inner_leaf
    
    upper_leaf = mav.Lattice(350., 350., 23., 2008, upper_leaf_comp, seed = 0)
    lower_leaf = mav.Lattice(350., 350., 0., 2008, lower_leaf_comp, seed = 0)
    upper_leaf.distribute()
    lower_leaf.distribute()
    bilayer = upper_leaf + lower_leaf
    
    vesicle + (bilayer.centroid() - vesicle.centroid())
    vesicle + np.array([0., 0., 300.])
    
    a2t = mav.Grid(240., 240., 0., 3, 3, a2t_comp, seed = 0)
    a2t.distribute() 
    a2t + (bilayer.centroid() - a2t.centroid())
    a2t + np.array([0., 0., 100.])
    
    vesi_bi_a2t = vesicle + bilayer + a2t
    
    vesi_bi_a2t.write_pdb("complexes/vesicle_bi_a2t.pdb")

'''
| Membrane-Vesicle Junction
'''
@generate_timed
def mem_vesi_benchmark():
    upper_leaf_comp   = mav.read_comp("compositions/upper_leaf_comp")
    lower_leaf_comp   = mav.read_comp("compositions/lower_leaf_comp")
    a2t_comp = mav.read_comp("compositions/a2t_comp")
    
    outer_leaf = mav.Sphere(125., 0., 3626, upper_leaf_comp, 
                        pore_radius = 20.0)
    inner_leaf = mav.Sphere(100., 0., 1857, lower_leaf_comp,
                        pore_radius = 20.0)
    
    outer_leaf.distribute()
    inner_leaf.distribute()
    
    vesicle  = outer_leaf + inner_leaf
    vesicle0 = copy.deepcopy(vesicle)
    vesicle0 + np.array([0., 0., 370.])
    vesi_vesi = vesicle + vesicle0

    a2t = mav.Grid(240., 240., 0., 3, 3, a2t_comp, seed = 0)
    a2t.distribute() 
    a2t + (vesi_vesi.centroid() - a2t.centroid())
    
    vesi_bi_a2t = vesi_vesi + a2t
    
    vesi_bi_a2t.write_pdb("complexes/vesi_vesi_a2t.pdb")

'''
| Covid Virion
'''

@generate_timed
def covid_benchmark():
    upper_leaf_comp   = mav.read_comp("compositions/upper_leaf_comp")
    lower_leaf_comp   = mav.read_comp("compositions/lower_leaf_comp")
    spike_comp = mav.read_comp("compositions/spike_comp")
    
    outer_leaf = mav.Sphere(500., 0., 51475, upper_leaf_comp, seed = 0)
    inner_leaf = mav.Sphere(475., 0., 46456, lower_leaf_comp, seed = 0)
    outer_leaf.distribute()
    inner_leaf.distribute()
    vesicle = outer_leaf + inner_leaf
    
    spike = mav.Sphere(600., 0., 75, spike_comp, seed = 0)
    spike.distribute()
    
    vesicle.remove_overlap(spike, 4.)
    
    covid = spike + vesicle
    
    covid.write_pdb("complexes/covid_viron.pdb")

if __name__ == "__main__":
    print("Bilayer")
    bilayer_benchmark()
    print("Porin Bilayer System")
    porin_gird_benchmark()
    print("Nanodisc")
    nanodisc_benchmark()
    print("Nanodisc Spike")
    nanodisc_spike_benchmark()
    print("Nanotube")
    nanotube_benchmark()
    print("Nanotube Array")
    nanotube_array_benchmark()
    print("Vesicle")
    vesicle_benchmark()
    print("Membrane-Vesicle Junction")
    mem_vesi_benchmark()
    print("Vesicle-Vesicle Junction")
    mem_vesi_benchmark()
    print("Covid")
    covid_benchmark()
