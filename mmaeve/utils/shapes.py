'''
Shapes
------

Code for generating the point distributions that are used for 
constructing systems.
'''

import numpy as np

def fib_lattice(length, width, height, n_points):
    ''' 
    Purpose:   Uses the fibonicci lattice to distribute points on the 
               surface of a unit square whose dimensions are then
               modified.
    Arguments: length) Float. Scale of lattice along the x-axis.
               width) Float. Scale of the lattice along the 
               y-axis.
               height) Float. z-coordinate of the points.
               n_points) Integer. Number of points to distribute.
    Returns:   N x 3 np.array of Floats. Points distributed on a
               lattice.
    '''
    if n_points == 1:
        return(np.array([[0., 0., height]]))
        
    eps = 0.5
    au_ratio = (1 + (5 ** 0.5)) / 2
    x = (np.arange(0, n_points, 1) / au_ratio) % 1 * length
    y = (np.arange(0, n_points, 1) + eps) \
        /                                 \
        ((n_points - 1) + (2. * eps)) * width
    z = np.repeat(0., n_points) + height
    
    return(np.stack([x, y, z], axis = -1))

def fib_disc(radius, height, n_points):
    ''' 
    Purpose:   Distribute points across the surface of a disc using 
               the fibonicci lattice.
    Arguments: radius) Float. Radius of the disc on which the points
               are distributed.
               height) Float. z-coordinate of the points.
               n_points) Integer. Number of points to distribute.
    Returns:   N x 3 np.array of Floats. Points distributed on a
               disc.
    '''
    lattice = fib_lattice(1., 1., height, n_points)
    theta = 2 * np.pi * lattice[:, 0]
    r = lattice[:, 1] ** 0.5

    x = r * radius * np.cos(theta)
    y = r * radius * np.sin(theta)
    z = lattice[:, 2]

    return(np.stack([x, y, z], axis = -1))

def fib_cylinder(radius, length, height, n_points):
    ''' 
    Purpose:   Distribute points across the surface of a cylinder 
               using the fibonicci lattice.
    Arguments: radius) Float. Radius of the cylinder.
               length) Float. Length of the cylinder.
               height) Float. Translation of the coordinates along the
               z-axis.
               n_points) Integer. Number of points to distribute.
    Returns:   N x 3 np.array of Floats. Points distributed on a
               cylinder.
    '''
    # Starting lattice
    lattice = fib_lattice(1., 1., 0., n_points)

    # Lattice => Cylindrical projection
    theta = 2 * np.pi * lattice[:, 0]
    y = lattice[:, 1]
    r = np.sqrt(y * y + 1.)

    # Cylindrical projection => Spherical coordinates
    theta = theta
    phi = np.arccos(y / r)
    r = r

    # Spherical coordinates => Cartesian coordinates
    x = r * np.sin(phi) * np.cos(theta) * radius
    y = r * np.sin(phi) * np.sin(theta) * radius
    z = r * np.cos(phi) * length

    return(np.stack([x, y, z], axis = -1))

def fib_sphere(radius, height, n_points):
    ''' 
    Purpose:   Distribute points across the surface of a sphere using 
               the fibonicc lattice.
    Arguments: radius) Float. The radius of the sphere on which points 
               are distributed.
               height) Float. Translation of the coordinates along the
               z-axis.
               num_points) Integer. The number of points to distribute 
               across the sphere's surface.
    Returns:   N x 3 np.array of Floats. Points distributed on the 
               surface of a sphere.
    '''
    # Starting lattice
    lattice = fib_lattice(1., 1., 0., n_points)

    # Lattice => Spherical projection
    theta = 2 * np.pi * lattice[:, 0]
    phi = np.arccos(1. - (2. * lattice[:, 1]))

    # Spherical projection => Cartesian coordinates
    x = np.sin(phi) * np.cos(theta) * radius
    y = np.sin(phi) * np.sin(theta) * radius
    z = np.cos(phi) * radius

    return(np.stack([x, y, z], axis = -1))

def circle(radius, height, n_points):
    '''
    Purpose:   Distribute points evenly on a circle.
    Arguments: num_points) Integer. Number of points to be 
               distributed.
               radius) Float. The desired radius of the circle.
               height) Float. Height of the system.
    Returns:   N x 3 np.array of Floats. Points distributed on a 
               circle.
    '''
    theta = np.linspace(0, 2 * np.pi, n_points, endpoint = False)
    x = radius * np.cos(theta)
    y = radius * np.sin(theta)
    z = np.repeat(height, len(x))
    
    return(np.stack([x, y, z], axis = -1))

def grid(length, width, height, l_number, w_number):
    '''
    Purpose:   Distribute l_number x w_number across a grid of 
               length x width dimensions.
    Arguments: length) Float. The desired size of the x-dimension.
               width) Float. The desired size of the y-dimension.
               height) Float. z-coordinate of the grid.
               l_number) Integer. Number of number to be distributed 
               along the x-axis.
               w_number) Integer. Number of number to be distributed 
               along the y-axis.
    Returns:   N x 3 np.array of Floats. Where each column 
               corresponds to the x, y, and z coordinates, 
               respectively.
    '''
    # Shift factor to center number in their respective grid 
    # subdivisions
    x_shift = (length / l_number) / 2.
    y_shift = (width / w_number) / 2.

    # Generate x and y coordinates
    x = np.arange(0, length, length / l_number) + x_shift
    y = np.arange(0, width, width / w_number) + y_shift
    
    # Create indicies to make combinations of x and y coordinates
    x_indicies = np.array([np.arange(0, l_number, 1)])
    x_indicies = np.repeat(x_indicies, w_number, axis = 0)
    x_indicies = x_indicies.flatten().astype(int)

    y_indicies = np.arange(0, w_number, 1)
    y_indicies = np.repeat(y_indicies, l_number, axis = 0).astype(int)
    
    # Concatenate and return points
    x = x[x_indicies]
    y = y[y_indicies]
    z = np.repeat(height, len(y_indicies))
    
    return(np.stack([x, y, z], axis = -1))
