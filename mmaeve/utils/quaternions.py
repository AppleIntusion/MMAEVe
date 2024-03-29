'''
Quaternions
-----------

Code related to quaternion operations.
'''

import numpy as np

def get_rot_quat(mobile, reference):
    '''
    Purpose:   Get the rotation and inverse rotation quaternions
               needed to align a given set of vectors with another 
               set of vectors.
    Arguments: mobile) N x 3 np.array of Floats. Each row is an
               individual vector. These points are considered
               mobile.
               reference) N x 3 np.array of Floats. Each row is an
               individual point. These points will be aligned to.
    Returns:   N x 3 np.array of Floats. Where each of the original 
               points has been aligned to a vector.
    '''
    # Check if vector norms are the same or inverse
    mobile_norm = mobile / \
                  np.linalg.norm(mobile, axis = 1)[:, np.newaxis] 
    ref_norm    = reference / \
                  np.linalg.norm(reference, axis = 1)[:, np.newaxis]
    dot = np.sum(mobile_norm * ref_norm, axis = 1)
    same_vec = dot >  0.9999
    inv_vec  = dot < -0.9999

    # Get the rotation axis and the angle to rotate about to align
    # "mobile" with "reference".
    # Numpy may throw a div by zero warning in the case of inverse
    # or aligned vectors. It will be handled.
    rot_axis = np.cross(mobile, reference)
    rot_axis = rot_axis / \
               np.reshape(np.linalg.norm(rot_axis, axis = 1), 
                          [len(rot_axis), 1]) # Normalize
    dot = np.sum(mobile * reference, axis = 1)
    angle = np.arccos(dot / \
              (np.linalg.norm(mobile, axis = 1)     * \
               np.linalg.norm(reference, axis = 1)
              )
            )

    # Handle the case of inverse vectors. A 180 degree rotation about 
    # the y-axis is performed.
    rot_axis[inv_vec, 0] = 0.
    rot_axis[inv_vec, 1] = 1.
    rot_axis[inv_vec, 2] = 0.
    angle[inv_vec] = np.pi

    # Create the rotation and inverse rotation quaternions.
    cos_half = np.cos(angle / 2.)
    sin_half = np.sin(angle / 2.)
    q0, q1, q2, q3 = cos_half,                  rot_axis[:, 0] * sin_half, \
                     rot_axis[:, 1] * sin_half, rot_axis[:, 2] * sin_half
    # Handle the case of aligned vectors. A zero rotation is performed.
    q0[same_vec] = 1.
    q1[same_vec] = 0.
    q2[same_vec] = 0.
    q3[same_vec] = 0.

    rot_quat = np.array([q0, q1, q2, q3])
    rot_quat = rot_quat.T
    inv_quat = np.array([q0, -q1, -q2, -q3])
    inv_quat = inv_quat.T

    return(rot_quat, inv_quat)

def quat_mult(quat0, quat1):
    '''
    Purpose:   Multiply two collections of quaternions.
    Arguments: quat0) N x 4 np.array of Floats. Each row is an
               individual quaternion.
               quat1) N x 4 np.array of Floats. Each row is an
               individual quaternion.
    Returns:   N x 4 np.array of Floats. Where each row 
               corresponds to a result of the individual quaternion 
               multiplications.
    '''
    w1, x1, y1, z1 = quat0[:, 0], quat0[:, 1], quat0[:, 2], quat0[:, 3]
    w2, x2, y2, z2 = quat1[:, 0], quat1[:, 1], quat1[:, 2], quat1[:, 3]
    return(
    np.column_stack([w1*w2 - x1*x2 - y1*y2 - z1*z2,
                     w1*x2 + x1*w2 + y1*z2 - z1*y2,
                     w1*y2 - x1*z2 + y1*w2 + z1*x2,
                     w1*z2 + x1*y2 - y1*x2 + z1*w2]))
