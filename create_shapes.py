import numpy as np

n=11
n2=11
radi = 1/(2*np.sin(np.pi/n))
radi2 = 1/(2*np.sin(np.pi/n2))

radii= [radi]
position = [0]
order = [n]
volume = [1]
width = [0.3]

for i, radius in enumerate(radii):
    for angle in range(order[i]):
        X = radius * np.sin(2*np.pi/order[i]*angle)
        Y = radius * np.cos(2*np.pi/order[i]*angle)
        Z = position[i]
        print('{0:10.5f} {1:10.5f} {2:10.5f}    {3:10.5f} {4:10.5f}'.format(X,Y,Z, volume[i], width[i]))