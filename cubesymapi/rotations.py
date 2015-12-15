import numpy as np

def rot_matrix_align(final, origin=(0,0,1)):
    origin = np.array(origin)/np.linalg.norm(np.array(origin))
    final = np.array(final)/np.linalg.norm(np.array(final))

    v = np.cross(origin, final)

    mat_v = np.array([[0, -v[2], v[1]],
                      [v[2], 0, -v[0]],
                      [-v[1], v[0], 0]])

    if np.linalg.norm(v) == 0:
        rot = np.identity(3)
    else:
        rot = np.identity(3) + mat_v + np.linalg.matrix_power(mat_v, 2)* (1-np.dot(origin, final))/(np.linalg.norm(v))
    return rot


def rot_matrix_z(angle):

    rot = np.array([[np.cos(angle), -np.sin(angle), 0],
                    [np.sin(angle), np.cos(angle), 0],
                    [0, 0, 1]])
    return rot


def rotate_align_z(coord, angle, center=(0, 0, 0), align=(0, 0, 1), radial=False):


    if radial:
        #Tranformation cilindrical -> cartessian
        x = coord[2] * np.cos(coord[1])
        y = coord[2] * np.sin(coord[1])
        z = coord[0] - center[2]

        coord = np.array([x, y, z])

        new_coord = np.dot(rot_matrix_align(align), np.dot(rot_matrix_z(angle), np.array(coord))) + np.array(center)


    else:
        new_coord =  np.dot(rot_matrix_align(align),
                      np.dot(rot_matrix_z(angle), np.array(coord) - np.array(center))
                      ) + np.array(center)



    return new_coord

