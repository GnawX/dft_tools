"""
for arbitrary many k points, calculate the ibz kpts
"""

from phonopy.interface.calculator import read_crystal_structure
import numpy as np
import spglib
import sys

# deprecated for loops
#from numba import njit
#@njit
#def ibzreduce(nk,nr,kf,kfr,tol=1.E-6):
#    kpt = []
#    wk = []
#    found = []
#    for i in range(nk):
#        if i not in found:
#            w = kf[i,3]
#            for j in range(i+1,nk):
#                if j not in found:
#                    lfound = False
#                    for k in range(nr):
#                        if max(abs(kf[i,0] - kfr[j,k,0]), \
#                               abs(kf[i,1] - kfr[j,k,1]), \
#                               abs(kf[i,2] - kfr[j,k,2])) < tol:
#                            lfound = True
#                            break
#                    if lfound:
#                        w += kf[j,3]
#                        found.append(j)
#            kpt.append(kf[i,0:3])
#            wk.append(w)
#    return kpt, wk

def get_rotationsk(atoms,time_reversal=True):
    """
    symmetry operations in fractional coords wrt k lattice
    """    
    rotations = spglib.get_symmetry(atoms)['rotations']
    nrotations = len(rotations)
    _rk = np.zeros((nrotations,3,3),dtype=np.int32)
    
    for i in range(nrotations):
        _rk[i,:,:] = np.linalg.inv(rotations[i].T)
    inv = np.eye(3,dtype=np.int32)*(-1)
    linv = False
    for i in range(nrotations):
        if np.array_equal(inv,rotations[i,:,:]):
            linv = True
            break
    if linv:
        rotationsk = _rk
    elif time_reversal:
        rotationsk = np.zeros((nrotations*2,3,3),dtype=np.int32)
        rotationsk[:nrotations,:,:] = _rk
        rotationsk[nrotations:,:,:] = np.tensordot(_rk,inv,(2,0))
    else:
        rotationsk = _rk

    return rotationsk


def writek(kfile,k,w):
    nk = len(k)
    with open(kfile,'w') as f:
         f.write('automatically kmesh \n')
         f.write('{:10d}\n'.format(nk))
         f.write('Reciprocal lattice \n')
         for i in range(nk):
             f.write(('{:20.14f}'*3 +'{:20.14f}\n').format(*k[i,:],w[i]))

def main():
   
    tol = 1.E-8 
    atoms, _ = read_crystal_structure('POSCAR')
    kpt = np.loadtxt('KPOINTS.DAT',skiprows=3)
    kf = kpt[:,0:3]
    wf = kpt[:,3]
    
    rotationsk = get_rotationsk(atoms)
    
    kfr = np.tensordot(rotationsk,kf,(2,1)).transpose(2,0,1)
    kfr = np.mod(kfr + 6.5 - 0.5*tol, [1.0,1.0,1.0]) - 0.5 + 0.5*tol
    kfr = np.around(kfr,decimals=8)
    nk,nr = np.shape(kfr)[0:2]

    ks = np.zeros((nk,3))

    for i in range(nk):
        ktmp = kfr[i,:,:]
        ks[i,:] = ktmp[np.lexsort(np.transpose(ktmp))][0]

    value, indices, inverse= np.unique(ks,axis=0,return_index=True,return_inverse=True)

    w = []
    for i in range(len(indices)):
        w.append(np.sum(wf[np.where(inverse==i)]))
    
    w = np.array(w)
    k = kf[indices]
    indk = np.lexsort(np.transpose(k))
    k = k[indk]
    w = w[indk]

    print(np.sum(w))
    writek('KPOINTS',k,w)

if __name__ == "__main__":
    main()
