from ase import atoms
from ase import io
import numpy as np

atoms = io.read('POSCAR')
cell = atoms.get_cell()
coor = atoms.get_scaled_positions()
#coor = np.loadtxt('POSCAR',skiprows=8)
csl = atoms.get_chemical_symbols()

strform = 3*'{:>20.16f}' + '\n'
strform2 = '{:<5s}' + 3*'{:>20.16f}' + '\n'
with open('qecoor','a') as f:
     f.write("CELL_PARAMETERS {angstrom}\n")
     for i in range(3):
         f.write(strform.format(*cell[i])) 
     f.write("ATOMIC_POSITIONS {crystal}\n")
     for i in range(len(csl)):
         f.write(strform2.format(csl[i],*coor[i]))
