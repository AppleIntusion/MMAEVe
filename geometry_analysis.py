import numpy as np

''' Purpose:   To calculate the smallest angle between two vectors.
    Arguments: (1) First vector.
               (2) Second vector.
               (3) Units which the angle value are returned in. The
                   default is 'rad' for radians and the other option is
                   'deg' for degrees.
'''
#def vec_angle(v1, v2, units='rad'):
#    np.seterr(all = 'raise')
#    #try:
#    #if (v1[0] == v1[1] == v1[2] == 0) or (v2[0] == v2[1] == v2[2] == 0):
#    #else:
#    angle = np.dot(np.array(v1), np.array(v2))
#    angle = angle / (np.linalg.norm(np.array(v1)) * np.linalg.norm(np.array(v2)))
#    angle = round(angle, 10)
#    angle = np.arccos(angle)
#        
#    #except FloatingPointError:
#    #angle = 0.0
#        
#    if units == 'rad':
#        pass
#    elif units == 'deg':
#        angle = (180 / np.pi)
#    return angle
#
#''' Purpose:   To calculate the angle between two vectors regardless of angle
#               size. You are measuring the angle between the first and the
#               second vectors. You must specify a normal vector that would result
#               if you traveled from v1 to v2. Be careful and make sure you
#               understand what it is actually doing if you want to use it. 
#               See S1.
#    Arguments: (1) First vector.
#               (2) Second vector.
#               (3) Vector normal to the first and second vectors.
#               (4) Units to return the angle in. Default is radians.
#'''
#def directional_angle(v1, v2, norm, units = "rad"):
#    angle = vec_angle(v1, v2)
#    trip_prod = np.dot(norm, np.cross(v1, v2))
#    if trip_prod < 0:
#        angle = (2 * np.pi) - angle
#    else:
#        pass
#    return angle
# 
#
#
#def rotate(v1, axis = 'x', theta = 0.):
#    ''' 
#    Purpose:   To rotate a vector by a specified amount around 
#               either the x, y, or z-axis.
#    Arguments: 1) Vector
#               2) String. The axis to rotate about: 'x', 'y', 
#                  or 'z'. 
#               3) Float. How far to rotate, in radians; default of 0.
#    Returns:   None. Updates Molecule instance.
#    '''
#    if axis == 'x':
#        v1[1], v1[2] =                       \
#                   (np.cos(theta) * v1[1]) - \
#                   (np.sin(theta) * v1[2])   \
#                   ,                         \
#                   (np.sin(theta) * v1[1]) + \
#                   (np.cos(theta) * v1[2])
#    if axis == 'y':
#        v1[0], v1[2] =                       \
#                   (np.cos(theta) * v1[0]) + \
#                   (np.sin(theta) * v1[2])   \
#                   ,                         \
#                   (np.cos(theta) * v1[2]) - \
#                   (np.sin(theta) * v1[0])   
#    if axis == 'z':
#        v1[0], v1[1] =                       \
#                   (np.cos(theta) * v1[0]) - \
#                   (np.sin(theta) * v1[1])   \
#                   ,                         \
#                   (np.sin(theta) * v1[0]) + \
#                   (np.cos(theta) * v1[1])   
def rotation_matrix_from_vectors(vec1, vec2):
    '''
    Purpose:   Find the rotation matrix that aligns vec1 to vec2
    Arguments: 1) vec1: A 3d "source" vector
               2) A 3d "destination" vector
    Returns:   A transform matrix (3x3) which when applied to vec1, 
               aligns it with vec2.
    '''
    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), \
           (vec2 / np.linalg.norm(vec2)).reshape(3)
    v = np.cross(a, b)
    if any(v): #if not all zeros then 
        c = np.dot(a, b)
        s = np.linalg.norm(v)
        kmat = np.array([[0, -v[2], v[1]], 
                            [v[2], 0, -v[0]], 
                            [-v[1], v[0], 0]])
        return(np.eye(3) + kmat + kmat.dot(kmat) * \
               ((1 - c) / (s ** 2)))
    else:
        # Cross of all zeros only occurs on identical directions
        return np.eye(3) 
