import numpy as np
import spglib
from ase import io
from scipy.io import FortranFile
import BoltzTraP2.sphere
import BoltzTraP2.fite
from ase.dft.kpoints import bandpath
import matplotlib.pyplot as plt


def read_eig(file):
    f = FortranFile(file,'r')
    f.read_record(dtype = np.int32)
    f.read_record(dtype = np.int32)
    f.read_record(dtype = np.int32)
    f.read_record(dtype = np.int32)
    f.read_record(dtype = np.int32)
    eig = f.read_record(dtype = np.float64)
    f.close()
    return eig.tolist()

mesh = [8,8,1]
atoms = io.read('POSCAR')

mapping, grid = spglib.get_ir_reciprocal_mesh(mesh,atoms)
kpoints = grid[np.unique(mapping)] / np.array(mesh,dtype=float)
ind = np.unique(mapping)+1
nk = len(ind)

f = FortranFile('q1/ACVKCAR','r')
nb = f.read_record(dtype=np.int32)[0]
f.close()

ebands = []
for i in range(len(ind)):
    f = 'q'+str(ind[i])+'/ACVKCAR'
    ebands += read_eig(f)
ebands = np.array(ebands).reshape((nb,nk),order='F')

class DFTdata:

      def __init__(self):
          self.kpoints = np.copy(kpoints)
          self.ebands = np.copy(ebands)
          self.mommat = None
       
      def get_lattvec(self):
          return atoms.get_cell().T

data = DFTdata()

lattvec = atoms.get_cell().T
#path = [(0,0,0),(0.5,0.0,0),(0.5,0.5,0),(0,0.5,0),(0,0,0)]
path = [(0,0,0),(0.5,0.0,0)]
bp = bandpath(path,atoms.cell,npoints=50)
kp = bp.kpts
(x, X, labels) = bp.get_linear_kpoint_axis()

equivalences = BoltzTraP2.sphere.get_equivalences(atoms,None,nk*8)
coeffs = BoltzTraP2.fite.fitde3D(data, equivalences)
eband, vband = BoltzTraP2.fite.getBands(kp,equivalences,lattvec,coeffs)

#np.savetxt('velocity',vband[0,:,20])

for i in range(50):
   plt.plot(x,eband[i,:])
plt.show()
