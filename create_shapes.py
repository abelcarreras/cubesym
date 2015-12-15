import numpy as np

radii= [4, 7]
position = [0, 0]
order = [3, 7]
volume = [2, 1]
width = [1.0, 1.3]

for i, radius in enumerate(radii):
    for angle in range(order[i]):
        X = radius * np.sin(2*np.pi/order[i]*angle)
        Y = radius * np.cos(2*np.pi/order[i]*angle)
        Z = position[i]
        print('{0:10.5f} {1:10.5f} {2:10.5f}    {3:10.5f} {4:10.5f}'.format(X,Y,Z, volume[i], width[i]))