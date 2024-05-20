import numpy as np
from math import sin, cos, pi, sqrt

n = 12 # number of atoms
c = 3 # lattice vector Ang

a = 1.42 # C-C

theta = 2*pi/n

#the radius
r = sqrt(a**2 - (c/n)**2)/2/sin(theta/2)

for i in range(n):
    x = r*cos(theta*i)
    y = r*sin(theta*i)
    z = c/n*i
    print(('{:15.6f}'*3).format(x,y,z))
