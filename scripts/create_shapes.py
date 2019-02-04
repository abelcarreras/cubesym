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


exit()
center = np.array([0.5, 0.5, 0.0])
vector = np.array([0, 0, 1])

cube = [[-1, -1, -1],
        [-1, -1, 1],
        [-1, 1, -1],
        [1, -1, -1],
        [1, 1, 1],
        [1, 1, -1],
        [1, -1, 1],
        [-1, 1, 1],
        ]

print("Distance         position")

for point in cube:
    d = np.linalg.norm(np.cross(vector, center-np.array(point)))/np.linalg.norm(vector)
 #   l = np.dot(vector, np.array(point))
    t = -np.dot(center-np.array(point), vector)/np.linalg.norm(vector)
    print("{}  {}".format(d, t))

exit()
