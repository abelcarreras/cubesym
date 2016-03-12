import numpy as np

radii= [1.15]
position = [0]
order = [5]
volume = [4]
width = [0.5]

for i, radius in enumerate(radii):
    for angle in range(order[i]):
        X = radius * np.sin(2*np.pi/order[i]*angle)
        Y = radius * np.cos(2*np.pi/order[i]*angle)
        Z = position[i]
        print('{0:10.5f} {1:10.5f} {2:10.5f}    {3:10.5f} {4:10.5f}'.format(X,Y,Z, volume[i], width[i]))