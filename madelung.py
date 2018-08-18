import numpy as np
from numpy.linalg import inv
from scipy.linalg import det
from math import exp, sqrt, erfc, pi
#calculate the madelung potential for orthorhombic lattice
#take into account the anisotropic dielectric tensor PRB 87, 094111(2013)

#parameters
a2b = 1.889725989 # ang to bohr
h2eV = 27.2114
sgm = 0.05 # unit of 1/Length
rm = 5
gm = 5
#real lattice vector in angstrom
r = np.array([31.6, 31.6, 31.6])
#
r = r*a2b
V = np.prod(r)
g = np.reciprocal(r)*2.0*pi
#dielectric tensor
eps = np.array([[3.4,0.0,0.0],[0.0,3.4,0.0],[0.0,0.0,3.4]])


#real sum
rsum = 0.0
for i in range(-rm,rm+1):
    for j in range(-rm,rm+1):
        for k in range(-rm,rm+1):
            ri = np.array([i,j,k])*r
            rii = np.dot(np.dot(ri,inv(eps)),ri)
            if rii > 0.0:
                rsum += 1/sqrt(det(eps))*erfc(sgm*sqrt(rii))/sqrt(rii)
#reciprocal sum
gsum = 0.0
for i in range(-gm,gm+1):
    for j in range(-gm,gm+1):
        for k in range(-gm,gm+1):
            gi = np.array([i,j,k])*g
            gii = np.dot(np.dot(gi,eps),gi)
            if gii > 0.0:
               gsum += exp(-gii/4/sgm**2)/gii*4*pi/V
#
vm = rsum + gsum - 2*sgm/sqrt(pi*det(eps)) - pi/V/sgm**2

print 'Madelung energy (eV):', vm/2*h2eV
