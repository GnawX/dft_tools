import numpy as np
from scipy.io import FortranFile

isp, ik, ib1, ib2 = input("Enter # of spin, kpt, cb, vb: ").split()
isp = int(isp); isp -= 1
ik = int(ik); ik -= 1
ib1 = int(ib1); ib1 -= 1
ib2 = int(ib2); ib2 -= 1

file = FortranFile('WAVEDER','r')
nb_tot, nbands_cder, nkpts, ispin = file.read_record(dtype= np.int32)
nodesn_i_dielectric_function = file.read_record(dtype= np.float64)
wplasmon = file.read_record(dtype= np.float64).reshape(3,3)
cder = file.read_record(dtype= np.complex64)
#cder = file.read_record(dtype= np.float32)
cder = np.reshape(cder, (nb_tot, nbands_cder, nkpts, ispin, 3), order='F')

print(cder[ib1,ib2,ik,isp,0])
print(cder[ib1,ib2,ik,isp,1])
print(cder[ib1,ib2,ik,isp,2])

#a = complex(0.)
#b = complex(0.)
#c = complex(0.)
#for i in range(864):
#   a = a + cder[i,i,0,0,0]
#   b = b + cder[i,i,0,0,1]
#   c = c + cder[i,i,0,0,2]
#print(a)
#print(b)
#print(c)
