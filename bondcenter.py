import numpy as np
from ase import io

#bond length cutoff
cut = 1.60  # A
strform = 3*'{:20.16f}'

atoms = io.read('POSCAR')
s = atoms.get_scaled_positions()

for i in range(len(atoms)):
    a = s[i]
    for j in range(i+1,len(atoms)):
        b = s[j]
        for k in range(3):
            if a[k]-b[k] > 0.5:
               b[k] = b[k] + 1
            elif a[k]-b[k] < -0.5:
               b[k] = b[k] - 1
        dcrys = a - b
        dcart = dcrys[0]*atoms.cell[0] + dcrys[1]*atoms.cell[1] + \
                dcrys[2]*atoms.cell[2]
        d = np.linalg.norm(dcart)
        if d < cut:
           mid = a/2 + b/2
           mid = np.mod(mid,1)
           print(strform.format(*mid))
