from scipy.interpolate import RegularGridInterpolator
from scipy import interpolate, optimize, integrate

#import multiprocessing as mp
import matplotlib.pyplot as plt
import numpy as np

import iofile
import rotations

def z_slides(x, y, z, function):
    return function([x, y, z])

def z_slides_r(y, x, z, function):
    return function([x, y, z])


def function_multi(angle, center, align, radial, ranges, function_use):
    x, y, z = ranges
    dens_rot = function_use(np.array([[[rotations.rotate_align_z([i, j, k], angle, center=center,
                                                                                   align=align,
                                                                                   radial=radial
                                                                 ) for k in z] for j in y] for i in x]))
    return {angle:dens_rot}

def function_measure(z_slide, ranges, function_use, epsabs, epsrel):
    minx, maxx, miny, maxy = ranges
    [integral, error] = integrate.dblquad(z_slides, miny, maxy, lambda x: minx, lambda x: maxx,
                                                    args=(z_slide, function_use),
                                                    epsabs=epsabs, epsrel=epsrel)

    return {z_slide: integral}


def function_measure_1d(z_coordinate, z_slide, ranges, function_use, epsabs, epsrel):
    minx, maxx, miny, maxy = ranges
    [integral, error] = integrate.quad(z_slides_r, miny, maxy,
                                                   args=(z_coordinate, z_slide, function_use),
                                                   epsabs=epsabs, epsrel=epsrel)

    return {z_slide: integral}


class Calculation:

    def __init__(self,
                 electronic_density,
                 ranges,
                 order=1,
                 center=(0, 0, 0),
                 align=(0, 0, 1),
                 radial=False,
                 r_grid=(20, 20, 20),
                 ):

        self._electronic_density = None
        self._order = order
        self._center = center
        self._align = align
        self._ranges_cart = ranges
        self._radial = radial

        self._density = None
        self._total_overlap = None
        self._measure = None

        self._grid = r_grid

        self._fn_electronic_density = None
        self._fn_overlap = None
        self._fn_density2 = None
        self._fn_density = None

        print('Input cube bounds (Cartesian)')
        print(' x: {0} {1}'.format(ranges[0][0], ranges[0][-1]))
        print(' y: {0} {1}'.format(ranges[1][0], ranges[1][-1]))
        print(' z: {0} {1}'.format(ranges[2][0], ranges[2][-1]))

        self._fn_electronic_density = RegularGridInterpolator(self._ranges_cart, electronic_density, bounds_error=False, fill_value=0)

        corner = [[ranges[0][0],  ranges[1][0],  ranges[2][0]],
                  [ranges[0][0],  ranges[1][0],  ranges[2][-1]],
                  [ranges[0][0],  ranges[1][-1], ranges[2][0]],
                  [ranges[0][-1], ranges[1][0],  ranges[2][0]],
                  [ranges[0][-1], ranges[1][-1], ranges[2][-1]],
                  [ranges[0][-1], ranges[1][-1], ranges[2][0]],
                  [ranges[0][-1], ranges[1][0],  ranges[2][-1]],
                  [ranges[0][0],  ranges[1][-1], ranges[2][-1]]]

        vector = self._align

        distances = [np.linalg.norm(np.cross(vector, center-np.array(point)))/np.linalg.norm(vector) for point in corner]
        positions = [-np.dot(center-np.array(point), vector)/np.linalg.norm(vector) for point in corner]
        max_rad = np.max(distances)

        print ('Oriented cube bounds')
        if radial:
            print('Using cylindrical')

    #        r1 = np.max([ranges[0][-1]-center[0],center[0]-ranges[0][1]])
    #        r2 = np.max([ranges[1][-1]-center[1],center[1]-ranges[1][1]])
    #        r3 = np.max([ranges[2][-1]-center[2],center[2]-ranges[2][1]])

        #    max_rad = np.max([r1, r2, r3])
        #    long_c = np.linspace(np.linalg.norm(center)-max_rad,np.linalg.norm(center)+max_rad, self._radial_grid[2])

            long_c = np.linspace(np.min(positions), np.max(positions), self._grid[2])
            angle = np.linspace(0, 2 * np.pi, self._grid[1])
            r = np.linspace(0, float(max_rad), self._grid[0])

            self._ranges = [long_c, angle, r]

            print(' z:     {0} {1} ({2})'.format(self._ranges[0][0], self._ranges[0][-1], self._grid[0]))
            print(' angle: {0} {1} ({2})'.format(self._ranges[1][0], self._ranges[1][-1], self._grid[1]))
            print(' r:     {0} {1} ({2})'.format(self._ranges[2][0], self._ranges[2][-1], self._grid[2]))

        else:
            print ('Using Cartesian')
            z = np.linspace(np.min(positions), np.max(positions), self._grid[2])
            y = np.linspace(center[1] - max_rad, center[1] + max_rad, self._grid[1])
            x = np.linspace(center[0] - max_rad, center[0] + max_rad, self._grid[0])

            self._ranges = [x, y, z]

            print(' x:     {0} {1} ({2})'.format(self._ranges[0][0], self._ranges[0][-1], self._grid[0]))
            print(' y:     {0} {1} ({2})'.format(self._ranges[1][0], self._ranges[1][-1], self._grid[1]))
            print(' z:     {0} {1} ({2})'.format(self._ranges[2][0], self._ranges[2][-1], self._grid[2]))

  #          self._ranges = self._ranges_cart

    def plot_full(self, step, rotation=0):
        x, y, z = self._ranges
        Y, X = np.meshgrid(y, x)
        density_slide = self.fn_electronic_density(np.array([[rotations.rotate_align_z([i, j, 0],
                                                                                       rotation,
                                                                                       center=self._center,
                                                                                       align=self._align
                                                                                       ) for j in y] for i in x]))

        for val in np.arange(self._ranges[0][0]+step, self._ranges[0][-1], step):
            density_slide += self.fn_electronic_density(np.array([[rotations.rotate_align_z([i, j, val],
                                                                                            rotation,
                                                                                            center=self._center,
                                                                                            align=self._align
                                                                                            ) for j in y] for i in x]))
        density_slide *=step
    #    levels = np.linspace(0, 0.5, 40)
        plt.xlabel("Bohr")
        plt.ylabel("Bohr")
        plt.contourf(X, Y, density_slide)
        plt.colorbar().set_label("Electronic density (e-/Bohr^2)")
        plt.show()


    def plot_full_rad(self, step, rotation=0):
        x, y, z = self._ranges
        Z, Y = np.meshgrid(z, y)
        density_slide = self.fn_electronic_density(np.array([[rotations.rotate_align_z([0, j, k],
                                                                                       rotation,
                                                                                       center=self._center,
                                                                                       align=self._align,
                                                                                       radial=self._radial
                                                                                       ) for k in z] for j in y]))

        for val in np.arange(self._ranges[0][0]+step, self._ranges[0][-1], step):
            density_slide += self.fn_electronic_density(np.array([[rotations.rotate_align_z([i, j, val],
                                                                                            rotation,
                                                                                            center=self._center,
                                                                                            align=self._align,
                                                                                            radial=self._radial
                                                                                            ) for j in y] for i in x]))
        density_slide *=step

        plt.contourf(Y, Z, density_slide)
        if self._radial:
            plt.xlabel("Angle (Radian)")
            plt.ylabel("Radius (Bohr)")
        else:
            plt.xlabel("Bohr")
            plt.ylabel("Bohr")
        plt.colorbar().set_label("Electronic density (e-/Bohr^3)")
        plt.show()

    def plot_slide(self, val, rotation=0):
        x, y, z = self._ranges

        Y, X = np.meshgrid(y, x)
        density_slide = self.fn_electronic_density(np.array([[rotations.rotate_align_z([i, j, val],
                                                                                       rotation,
                                                                                       center=self._center,
                                                                                       align=self._align
                                                                                       ) for j in y] for i in x]))
    #    levels = np.linspace(0, 0.5, 40)
        plt.xlabel("Bohr")
        plt.ylabel("Bohr")
        plt.contourf(X, Y, density_slide)
        plt.colorbar().set_label("Electronic density (e-/Bohr^3)")
        plt.show()

    def plot_slide_rad(self, val, rotation=0):
        x, y, z = self._ranges
        Z, Y = np.meshgrid(z, y)
        density_slide = self.fn_electronic_density(np.array([[rotations.rotate_align_z([val, j, k],
                                                                                        rotation,
                                                                                        center=self._center,
                                                                                        align=self._align,
                                                                                        radial=self._radial
                                                                                        ) for k in z] for j in y]))

        plt.contourf(Y, Z, density_slide)
        if self._radial:
            plt.xlabel("Angle (Radian)")
            plt.ylabel("Radius (Bohr)")
        else:
            plt.xlabel("Bohr")
            plt.ylabel("Bohr")
        plt.colorbar().set_label("Electronic density (e-/Bohr^3)")
        plt.show()


    def get_density(self):
        if self._density is None:
            x, y, z = self._ranges
            self._density = self.fn_electronic_density(np.array([[[rotations.rotate_align_z([i, j, k],
                                                                                  0,
                                                                                  center=self._center,
                                                                                  align=self._align,
                                                                                  radial=self._radial,
                                                                                  ) for k in z] for j in y] for i in x]))
        return self._density

    def get_total_overlap(self):
        import multiprocessing

        if self._total_overlap is None:

            overlappings = {}

            def log_result(result):
                overlappings.update(result)

            pool = multiprocessing.Pool(processes=max(multiprocessing.cpu_count()-1, 1))
     #       print('using {} processes'.format(pool._processes))
            for i in range(self._order):
                angle = 2*np.pi/self._order * i
                pool.apply_async(function_multi,
                                 args=(angle, self._center, self._align, self._radial, self._ranges, self.fn_electronic_density,),
                                 callback=log_result)
            pool.close()
            pool.join()

            overlappings = [np.multiply(self.get_density(),overlappings[k]) for k in sorted(overlappings)]

            self._total_overlap = sum(overlappings)/self._order

        return self._total_overlap

    def get_total_measure(self,  epsabs=1e2, epsrel=1e2):
# Not working yet
        ranges = [self._ranges[0][0],
                  self._ranges[0][-1],
                  self._ranges[1][0],
                  self._ranges[1][-1],
                  self._ranges[2][0],
                  self._ranges[2][-1]]

        minx, maxx, miny, maxy, minz, maxz = ranges
        print ranges

        print z_slides(0,0,0,self.fn_density)

        area = integrate.tplquad(z_slides, minx, maxx, lambda y: miny, lambda y: maxy, lambda y, z: minz, lambda y, z: maxz,
                                 args=(self.fn_density),
                                 epsabs=epsabs, epsrel=epsrel)
        print area

    def get_measure(self, n_points=None, epsabs=1e2, epsrel=1e2, measure_error=1E-5):
        import multiprocessing

        if self._measure is None:

            if n_points:
                measure_points = np.linspace(self._ranges[2][0], self._ranges[2][-1], n_points)
            else:
                measure_points = np.linspace(self._ranges[2][0], self._ranges[2][-1], self._grid[2])

            ranges = [self._ranges[0][0],
                      self._ranges[0][-1],
                      self._ranges[1][0],
                      self._ranges[1][-1]]

            print('#     coordinate     measure(C{0})       overlap        density2        density'.format(self._order))
            measure = {'symmetry' : [], 'coordinate' : [], 'overlap' : [], 'density2' : [], 'density' : []}

            overlap = {}
            density = {}
            density2 = {}

            def log_overlap(result):
                overlap.update(result)

            def log_density(result):
                density.update(result)

            def log_density2(result):
                density2.update(result)

            pool = multiprocessing.Pool(processes=max(multiprocessing.cpu_count()-1, 1))
        #       print('using {} processes'.format(pool._processes))
            for z_slide in measure_points:
                pool.apply_async(function_measure,
                                 args=(z_slide, ranges, self.fn_overlap, epsabs, epsrel),
                                 callback=log_overlap)

                pool.apply_async(function_measure,
                                 args=(z_slide, ranges, self.fn_density2, epsabs, epsrel),
                                 callback=log_density2)

                pool.apply_async(function_measure,
                                 args=(z_slide, ranges, self.fn_density, epsabs, epsrel),
                                 callback=log_density)

            pool.close()
            pool.join()

            overlap = [overlap[k] for k in sorted(overlap)]
            density2 = [density2[k] for k in sorted(density2)]
            density = [density[k] for k in sorted(density)]

            for integral_overlap, integral_density2, integral_density, z_slide in zip(overlap, density2, density, measure_points):
                if abs(integral_density2-integral_overlap) < measure_error or integral_density2 == 0:
                    normalized_measure = 0
                else:
                    normalized_measure = (1-integral_overlap/integral_density2)*100

                if self._radial:
                    integral_density *= z_slide
                    integral_density2 *= z_slide
                    integral_overlap *=z_slide

                print('{0:15.8f} {1:15.8f} {2:15.8f} {3:15.8f} {4:15.8f}'.format(z_slide,
                                                                                 normalized_measure,
                                                                                 integral_overlap,
                                                                                 integral_density2,
                                                                                 integral_density))

                measure['symmetry'].append(normalized_measure)
                measure['coordinate'].append(z_slide)
                measure['density2'].append(integral_density2)
                measure['density'].append(integral_density)
                measure['overlap'].append(integral_overlap)


            # Total measure
            total_overlap = integrate.simps(measure['overlap'], measure['coordinate'],even='avg')
            total_density2 = integrate.simps(measure['density2'], measure['coordinate'],even='avg')
            total_density = integrate.simps(measure['density'], measure['coordinate'],even='avg')
            total_symmetry = (1-total_overlap/total_density2)*100
            measure.update({'total_symmetry' : total_symmetry})
            print ('\nTotal measure: {0}'.format(total_symmetry))
            print ('\nTotal density: {0}'.format(total_density))

            self._measure = measure

        return self._measure

    def get_measure_1D(self, n_points=None, epsabs=1e2, epsrel=1e2, measure_error=1E-5, z_coordinate=0):

        if self._measure is None:

            if n_points:
                measure_points = np.linspace(self._ranges[2][0], self._ranges[2][-1], n_points)
            else:
                measure_points = np.linspace(self._ranges[2][0], self._ranges[2][-1], self._grid[2])

  #          minx = self._ranges[0][0]
  #          maxx = self._ranges[0][-1]

            miny= self._ranges[1][0]
            maxy = self._ranges[1][-1]
       #     print(minx, maxx)
       #     print(miny, maxy)

            print('#     coordinate     measure(C{0})       overlap        density2        density'.format(self._order))
            measure = {'symmetry' : [], 'coordinate' : [], 'overlap' : [], 'density2' : [], 'density' : []}
            for z_slide in measure_points:
                [integral_overlap, error_overlap] = integrate.quad(z_slides_r,
                                                                   miny, maxy,
                                                                   args=(z_coordinate, z_slide , self.fn_overlap),
                                                                   epsabs=epsabs, epsrel=epsrel)

                [integral_density2, error_density2] = integrate.quad(z_slides_r,
                                                                     miny, maxy,
                                                                     args=(z_coordinate, z_slide, self.fn_density2),
                                                                     epsabs=epsabs, epsrel=epsrel)

                [integral_density, error_density] = integrate.quad(z_slides_r,
                                                                   miny, maxy,
                                                                   args=(z_coordinate, z_slide, self.fn_density),
                                                                   epsabs=epsabs, epsrel=epsrel)

                if abs(integral_density2-integral_overlap) < measure_error:
                    normalized_measure = 0
                else:
                    try:
                        normalized_measure = (1-integral_overlap/integral_density2)*100
                    except ZeroDivisionError:
                        normalized_measure = 0.0

                if self._radial:
                    integral_density *= z_slide
                    integral_density2 *= z_slide
                    integral_overlap *= z_slide

                print('{0:15.8f} {1:15.8f} {2:15.8f} {3:15.8f} {4:15.8f}'.format(z_slide,
                                                                                 normalized_measure,
                                                                                 integral_overlap,
                                                                                 integral_density2,
                                                                                 integral_density))

                measure['symmetry'].append(normalized_measure)
                measure['coordinate'].append(z_slide)
                measure['density2'].append(integral_density2)
                measure['density'].append(integral_density)
                measure['overlap'].append(integral_overlap)

            # Total measure
            total_overlap = integrate.simps(measure['overlap'], measure['coordinate'],even='avg')
            total_density2 = integrate.simps(measure['density2'], measure['coordinate'],even='avg')
            total_density = integrate.simps(measure['density'], measure['coordinate'],even='avg')
            total_symmetry = (1-total_overlap/total_density2)*100
            measure.update({'total_symmetry' : total_symmetry})
            print ('\nTotal measure: {0}'.format(total_symmetry))
            print ('\nTotal density: {0}'.format(total_density))

            self._measure = measure

        return self._measure

    def plot_measure(self):

        fig, ax1 = plt.subplots()

        ax2 = ax1.twinx()
        ax3 = ax1.twiny()

        ax1.set_xlabel('Bohr')
        ax1.set_ylabel('Asymmetry %')
        ax2.set_ylabel('Density (e-/Bohr)')
        ax3.set_xlabel('Electrons')

        ax3.xaxis.grid()

        data = self.get_measure()

        step = data['coordinate'][1] - data['coordinate'][0]

        fn = interpolate.interp1d(data['coordinate'], np.cumsum(data['density'])*step ,bounds_error=False)
        def fun(x, n):
            return fn(x) - n

        total_electron = int(fun(data['coordinate'][-1], 0))

        n_step = abs(total_electron/10)+1
        positions = [optimize.bisect(fun, data['coordinate'][0], data['coordinate'][-1], args=(i, )) for i in range(1, total_electron+1, n_step)]

        positions = (np.array(positions) - data['coordinate'][0])/(data['coordinate'][-1] - data['coordinate'][0])
        plt.xticks(positions, range(1, total_electron+1, n_step))

        asim_p, = ax1.plot(data['coordinate'], data['symmetry'],
                           label='Asymmetry C{0}'.format(self._order), color='red', linewidth=4)
        #ax2.plot(data['coordinate'], data['density2'], label='Density2')
        dens_p, =ax2.plot(data['coordinate'], data['density'], label='Density', linestyle='--')
        #ax2.plot(data['coordinate'], data['overlap'], label='Overlap')
        plt.ylim(0,100)
        plt.legend(handles=[asim_p, dens_p],loc=1)

        plt.show()

    #Properties
    @property
    def fn_electronic_density(self):
        return self._fn_electronic_density

    @property
    def fn_overlap(self):
        if self._fn_overlap is None:
            x, y, z = self._ranges
            total_overlap = self.get_total_overlap()
            self._fn_overlap = RegularGridInterpolator((x, y, z), total_overlap, bounds_error=False)
        return self._fn_overlap

    @property
    def fn_density2(self):
        if self._fn_density2 is None:
            x, y, z = self._ranges

            density = self.get_density()
            self._fn_density2 = RegularGridInterpolator((x, y, z), np.square(density), bounds_error=False)
        return self._fn_density2

    @property
    def fn_density(self):
        if self._fn_density is None:
            x, y, z = self._ranges

            density = self.get_density()
            self._fn_density = RegularGridInterpolator((x, y, z), density, bounds_error=False)
        return self._fn_density


if __name__ == "__main__":
    #ranges, cube_density = get_density_cube('C3B9.cube')

    ranges, cube_density = iofile.get_density_gaussian('example/potential')
    #ranges, cube_density = get_density_gaussian()

    calculation = Calculation(cube_density, ranges, order=3, align=[0, 0, 1], center=[0, 0, 0], radial=True)
    #calculation = Calculation(cube_density, ranges)

    data = []
    for z in np.arange(-10, 10, 0.1):
        print (z, calculation.fn_density([z, 0, 0]))
        data.append(calculation.fn_density([0, 0, z])[0])
    #plt.plot(np.arange(-10, 10, 0.1), data)

    calculation.plot_slide(2, rotation=0)

    calculation.plot_slide_rad(0, rotation=0)

    calculation.plot_measure()
