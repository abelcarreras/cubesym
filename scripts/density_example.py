from cubesymapi import iofile, Calculation
import numpy as np
import matplotlib.pyplot as plt


# read gaussian cube file
ranges, cube_density = iofile.get_density_cube('../data/C3B9.cube')

# create Calculation instance
calculation = Calculation(cube_density, ranges)

# get density from [x:0, y:0, z:-6] to [x:0, y:0, z:6]
density_list = []
for z in np.arange(-6, 6, 0.1):
    den = calculation.fn_electronic_density([0, 0, z])[0]
    print ('point: {:.2f} {:.2f} {:.2f}'.format(*[0, 0, z]) + ' density: {}'.format(den))
    density_list.append(den)

# plot density
plt.plot(np.arange(-6, 6, 0.1), density_list)
plt.xlabel('Z direction [Bohr]')
plt.ylabel('electronic density [$Bohr^{-3}$]')
plt.show()

# plot slide at z=4.0 (for reference)
calculation.plot_slide(4.0, rotation=0)
