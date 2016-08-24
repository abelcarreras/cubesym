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


#A little faster
def get_density_gaussian(file):
    gaussian_file = open(file, "r")

    gaussian_data = []
    for line in gaussian_file.readlines():
        if line.find('#') == 0:
            continue
        gaussian_data.append(line.split())

    limit_data = np.array(gaussian_data[:3], dtype=float)
    x = linspace(*limit_data[0, :])
    y = linspace(*limit_data[1, :])
    z = linspace(*limit_data[2, :])
    V = zeros((limit_data[0, 2], limit_data[1, 2], limit_data[2, 2]))
    V_temp = zeros((limit_data[0, 2], limit_data[1, 2], limit_data[2, 2]))

    potential_data = np.array(gaussian_data[3:], dtype=float)

    for potential in potential_data:
        V_temp.fill(0)
        for i, cx in enumerate(x):
            for j, cy in enumerate(y):
                for k, cz in enumerate(z):
                    V_temp[i, j, k] = gauss3d([cx, cy, cz], potential[:3], potential[3], potential[4])
        V += V_temp
    return [x, y, z], V

def get_density_gaussian_old(file):
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
    import mmap

    with open (file, "r+") as f:
        cube_file = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
    #    position_number=file_map.find('NIONS =')
    #    file_map.seek(position_number+7)
        cube_file.readline()
        cube_file.readline()
        line = cube_file.readline()

        natom = int(line.split()[0])
        xmin, ymin, zmin  = [float(string) for string in line.split()][1:4]

        line = cube_file.readline()
        nx = int(line.split()[0])
        xstep = float(line.split()[1])
        line = cube_file.readline()
        ny = int(line.split()[0])
        ystep = float(line.split()[2])
        line = cube_file.readline()
        nz = int(line.split()[0])
        zstep = float(line.split()[3])

        for i in range(abs(natom)):
            cube_file.readline()

        V = np.zeros([nx, ny, nz])

        if int(natom) < 0:

            line = cube_file.readline()
            n_orbitals = int(line.split()[0])

            nl = int(np.ceil((n_orbitals+1)/10.))
            orbitals = line.split()[1:]

            for i in range(1,nl):
                line = cube_file.readline()
                orbitals += (line.split())

            print("Reading MO: {0}".format(' '.join(orbitals)))
            w = 2

        else:
            w = 1
            n_orbitals = 1
            print("Reading Density")

        for i in range(nx):

            line = cube_file.readline()
            density = line.split()
            while len(density) < ny*nz*n_orbitals:
                line = cube_file.readline()
                density += line.split()

            for j in range(ny):
                for k in range(nz):
                    for mo in range(n_orbitals):
                        V[i, j, k] += pow(float(density[j*nz*n_orbitals + k*n_orbitals + mo]), w)

    x = linspace(xmin, xmin+xstep*nx, nx)
    y = linspace(ymin, ymin+ystep*ny, ny)
    z = linspace(zmin, zmin+zstep*nz, nz)

    return [x, y, z], V


def get_density_cube_old(file):
    cube_file = open(file, "r")
    cube_data = cube_file.readlines()
    natom = np.abs(int(cube_data[2].split()[0]))
    xmin, ymin, zmin = np.array(cube_data[2].split()[1:4],dtype=float)
    nx, xstep = int(cube_data[3].split()[0]), float(cube_data[3].split()[1])
    ny, ystep = int(cube_data[4].split()[0]), float(cube_data[4].split()[2])
    nz, zstep = int(cube_data[5].split()[0]), float(cube_data[5].split()[3])

    V = np.zeros([nx, ny, nz])

    print(natom)
    if int(cube_data[2].split()[0]) < 0:
        n_orbitals = int(cube_data[6+natom].split()[0])
        nl = int(np.ceil((n_orbitals+1)/10.))
        orbitals = cube_data[6+natom].split()[1:]
        for i in range(1,nl):
            orbitals += (cube_data[i+6+natom].split())

        print("Reading MO: {0}".format(' '.join(orbitals)))
        density_data = " ".join(cube_data[nl+6+natom:]).split()
        w = 2

        for mo in range(n_orbitals):
            for i in range(nx):
                for j in range(ny):
                    for k in range(nz):
                        V[i, j, k] += pow(float(density_data[i*ny*nz*n_orbitals + j*nz*n_orbitals + k*n_orbitals + mo]), w)


    else:
        print("Reading Density")
        density_data = " ".join(cube_data[6+natom:]).split()
        w = 1

        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    V[i, j, k] = pow(float(density_data[i*ny*nz + j*nz +k]), w)

    cube_file.close()
    x = linspace(xmin, xmin+xstep*nx, nx)
    y = linspace(ymin, ymin+ystep*ny, ny)
    z = linspace(zmin, zmin+zstep*nz, nz)

    return [x, y, z], V

if __name__ == "__main__":

    get_density_cube('/Users/abel/CIRCULENS/nous_anells/c5h5_pi.cube')
  #  get_density_cube('/Users/abel/CIRCULENS/nous_anells/c5h5.cube')