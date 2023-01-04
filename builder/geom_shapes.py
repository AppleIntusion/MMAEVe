'''-------------------------------------------------------*
| Title:   Geometric Shapes                               |
|                                                         |
| Author:  Gubbin Eel (Satanic Overlord of the Swamp)     |
|                                                         |
| Purpose: Defines functions used to assign points to one |
|          of three geometric scaffolds: a sphere, nested |
|          sphere, or two rectangles.                     |
*-------------------------------------------------------'''

import numpy as np

''' Purpose:   To calculate the smallest angle between two vectors.
    Arguments: (1) First vector.
               (2) Second vector.
               (3) Units which the angle value are returned in. The
                   default is 'rad' for radians and the other option is
                   'deg' for degrees.
'''
def vec_angle(v1, v2, units='rad'):
    np.seterr(all = 'raise')
    #try:
    #if (v1[0] == v1[1] == v1[2] == 0) or (v2[0] == v2[1] == v2[2] == 0):
    #else:
    angle = np.dot(np.array(v1), np.array(v2))
    angle = angle / (np.linalg.norm(np.array(v1)) * np.linalg.norm(np.array(v2)))
    angle = round(angle, 10)
    angle = np.arccos(angle)
        
    #except FloatingPointError:
    #angle = 0.0
        
    if units == 'rad':
        pass
    elif units == 'deg':
        angle = (180 / np.pi)
    return angle

''' Purpose:   To calculate the angle between two vectors regardless of angle
               size. You are measuring the angle between the first and the
               second vectors. You must specify a normal vector that would result
               if you traveled from v1 to v2. Be careful and make sure you
               understand what it is actually doing if you want to use it. 
               See S1.
    Arguments: (1) First vector.
               (2) Second vector.
               (3) Vector normal to the first and second vectors.
               (4) Units to return the angle in. Default is radians.
'''
def directional_angle(v1, v2, norm, units = "rad"):
    angle = vec_angle(v1, v2)
    trip_prod = np.dot(norm, np.cross(v1, v2))
    if trip_prod < 0:
        angle = (2 * np.pi) - angle
    else:
        pass
    return angle
 


''' Purpose:   To evenly distribute points across a sphere with a
               specified radius are return a list of coordinates
               corresponding to those points.
    Arguments: (1) The number of points to distribute across the sphere's
                   surface.
               (2) The radius of the sphere which points should be
                   distributed on.
'''

def fib_sphere(num_points, radius = 1.0):
    ga = (3 - np.sqrt(5)) * np.pi # Value of the golden angle.
    
    # Golden angle increments along the number of points.
    ang = ga * np.arange(num_points)

    # Create a range of points from -1 to 1 to create a unit circle.
    z = np.linspace((1 / num_points) - 1, 1 - (1 / num_points), num_points)

    # A list of the radii at each height step of the unit circle.
    rad = np.sqrt(1 - (z * z))

    # Determine x and y values given angles.
    x = rad * np.cos(ang)
    y = rad * np.sin(ang)

    x = np.transpose([list(x)])
    y = np.transpose([list(y)])
    z = np.transpose([list(z)])
    xyz = np.concatenate([x, y, z], axis = 1)
    xyz = xyz * radius
    #x = np.array(x) * radius 
    #y = np.array(y) * radius 
    #z = np.array(z) * radius 
    rad = np.sqrt((x * x) + (y * y) + (z * z))

    return np.array(xyz)

def sun_radius(k, num_points, boundary_points):
    if k > num_points - boundary_points:
        r = 1
    else:
        r = np.sqrt(k - 0.5) / np.sqrt(num_points - (boundary_points + 1) / 2)
    return r

def sunflower(num_points, radius, alpha = 2.0, height = 0.0):
    '''
    Purpose   :     This function uses the sunflower distribution to evenly distribute
            points within a circle of a given radius
    Arguments :    1) Num_points --> points to be distributed
            2) Alpha --> important for populating. 2 is smoother than 0
            3) Radius --> use your brain dipshit
    Returns   : List of lists of x,y coords. [[x0, y0]...[xn, yn]]
    '''
    boundary_points = round(alpha * np.sqrt(num_points))
    phi = (np.sqrt(5.0) + 1.0) * 0.5
    point_list = []
    for k in list(range(1, num_points + 1)):
        
        r = sun_radius(k, num_points, boundary_points) * radius
        theta = 2.0 * np.pi * k / phi ** 2.0
        x, y = r * np.cos(theta), r * np.sin(theta)
        point_list.append([x, y, height])
    return point_list

def circle(num_points, radius, height):
    point_list = []
    theta = np.linspace(0, 2 * np.pi, num_points, endpoint=False)
    
    for angle in theta:
        x = radius * np.cos(angle) 
        y = radius * np.sin(angle)
        point_list.append([x, y, height])
    return point_list

def grid(num_points, length, width, height = 0.0):
    point_collection = []

    if num_points == 0:
        return point_collection

    n = num_points
    l = length
    w = width  
    if width != 0. or length != 0.:
        yn = np.sqrt(((w * n) / l) +                      \
                     (((w - l) ** 2) / (4 * (l ** 2)))) + \
             ((w - l) / (2 * l))
    else:
        yn = 1
    
    yn = int(yn)
    xn = int(num_points / yn)

    xmin = 0. - (length / 2.)
    xmax = (length / 2.)
    ymin = 0. - (width / 2.)
    ymax = (width / 2.)

    x_coords = np.linspace(xmin, xmax, xn)
    y_coords = np.linspace(ymin, ymax, yn)
    for ii in range(len(x_coords)):
        for jj in range(len(y_coords)):
            x = x_coords[ii]
            y = y_coords[jj]
            point_collection.append([x, y, height])
    point_collection = np.array(point_collection)
    return point_collection

'''--------------
| Supplementary |
--------------'''

'''
(S1) If you just solve for an angle using the cosine you get the smaller
     of the two possible angles. This method allows you to specify which
     of the two angles you want to use by taking a triple product.
     You specify two angles v1 and v2. This method will give you the
     angle from v1 to v2 even if it is larger than 180 degrees. You do
     this by specifying "norm" which is a vector normal to v1 and v2 and
     would result from their cross product. Then the cross product of v1
     and v2 is dotted with this normal. If it is negative then you know
     that you have an angle greater than 180 degrees and that will be
     handled. CAUTION: THIS FUNCTION WILL ONLY PROVIDE THE ANGLE IF YOU
     USE IT IN THE WAY SPECIFIED.
'''
