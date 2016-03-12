from numpy import linspace, zeros, array
import numpy as np

import warnings
warnings.filterwarnings('ignore')


def gauss3d(coor, center, volume, width):
    a = volume/(width * np.sqrt(2*np.pi))**3

    return a * (np.exp(-(pow(coor[0]-center[0],2)/(2*pow(width, 2)) +
                         pow(coor[1]-center[1],2)/(2*pow(width, 2)) +
                         pow(coor[2]-center[2],2)/(2*pow(width, 2)))))



def get_density_radial_test():
    x = linspace(-15, 15, 60)
    y = linspace(-15, 15, 60)
    z = linspace(-15, 15, 60)
    V = zeros((60, 60, 60))

    rot = 3
    width = 2.4

    for i, cx in enumerate(x):
        for j, cy in enumerate(y):
            for k, cz in enumerate(z):
                phi = np.arctan2(cy, cx)
                rad = np.sqrt(pow(cx,2)+pow(cy,2))
                V[i, j, k] = pow(np.sin(phi*rot/2),2) * np.exp(-(pow(cz,2)/(2*pow(width, 2)))) * np.exp(-(pow(rad,2)/(2*pow(width*1.5, 2))))

    return [x, y, z], V



def get_density_gaussian(file):
    gaussian_file = open(file, "r")

    gaussian_data = []
    for line in gaussian_file.readlines():
        if line.find('#') == 0:
            continue
        gaussian_data.append(line.split())

    limit_data = np.array(gaussian_data[:3] , dtype=float)
    x = linspace(*limit_data[0, :])
    y = linspace(*limit_data[1, :])
    z = linspace(*limit_data[2, :])
    V = zeros((limit_data[0, 2], limit_data[1, 2], limit_data[2, 2]))

    potential_data = np.array(gaussian_data[3:], dtype=float)

    for i, cx in enumerate(x):
        for j, cy in enumerate(y):
            for k, cz in enumerate(z):
                for potential in potential_data:
                    V[i, j, k] += gauss3d([cx, cy, cz], potential[:3], potential[3], potential[4])

    return [x, y, z], V


def get_density_cube(file):
    cube_file = open(file, "r")
    cube_data = cube_file.readlines()
    natom = np.abs(int(cube_data[2].split()[0]))
    xmin, ymin, zmin = np.array(cube_data[2].split()[1:4],dtype=float)
    nx, xstep = int(cube_data[3].split()[0]), float(cube_data[3].split()[1])
    ny, ystep = int(cube_data[4].split()[0]), float(cube_data[4].split()[2])
    nz, zstep = int(cube_data[5].split()[0]), float(cube_data[5].split()[3])

    print(natom)
    if int(cube_data[2].split()[0]) < 0:
        print("Reading MO {0} {1}".format(*cube_data[6+natom].split()))
        density_data = " ".join(cube_data[6+natom+1:]).split()
        w = 2

    else:
        print("Reading Density")
        density_data = " ".join(cube_data[6+natom:]).split()
        w = 1

    V = np.zeros([nx, ny, nz])

    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                V[i, j, k] = pow(float(density_data[i*ny*nz + j*nz +k]), w)

    cube_file.close()
    x = linspace(xmin, xmin+xstep*nx, nx)
    y = linspace(ymin, ymin+ystep*ny, ny)
    z = linspace(zmin, zmin+zstep*nz, nz)

    return [x, y, z], V


