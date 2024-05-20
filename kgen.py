from phonopy.interface.calculator import read_crystal_structure
import numpy as np
import spglib


mesh  = [4, 4, 1]
shift = [0, 0, 0]

atoms, _ = read_crystal_structure('POSCAR')
mapping, grid = spglib.get_ir_reciprocal_mesh(mesh, atoms, shift)

kpt, counts = np.unique(mapping,return_counts=True)
nkpts = len(kpt)
kpt = grid[kpt] / np.array(mesh, dtype=float)

with open('KPOINTS.DAT','w') as f:
     f.write('automatically kmesh \n')
     f.write('{:10d}\n'.format(nkpts))
     f.write('Reciprocal lattice \n')
     for i in range(nkpts):
         f.write(('{:20.14f}'*3 +'{:13d}\n').format(*kpt[i,:],counts[i]))
