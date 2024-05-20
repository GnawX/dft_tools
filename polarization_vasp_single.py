import numpy as np
from pymatgen.io.vasp.outputs import Outcar
from pymatgen import Structure
from pymatgen.analysis.ferroelectricity.polarization import get_total_ionic_dipole

e    =  1602.176634 # e/A2 --> \muC/cm2

# modulo operator
def mod(a,b):
    f = np.fmod(a,b)
    for i in range(len(f)):
        if abs(abs(f[i]) - b[i]) < 1.0e-5:
           f[i] = 0.0
    return f

# file readin
o = Outcar('OUTCAR')
s = Structure.from_file('POSCAR')
l = s.lattice

pe = o.p_elec # cartesian
pe = l.get_vector_along_lattice_directions(pe) # crystal
# ionic part from vasp
pi = o.p_ion # cartesian
pi = l.get_vector_along_lattice_directions(pi) # crystal
# ionic part from pymatgen
pi = get_total_ionic_dipole(s,o.zval_dict)
pi = -pi
# polarization quantum
pq = np.array(l.abc)
# modulo
pe = mod(pe,pq)
pi = mod(pi,pq)
p = pe + pi
p = mod(p,pq)
# convert e/A2 to muC/cm2
p = p/l.volume*e
pq = pq/l.volume*e
# calculate the norm
pcart = 0.0
for i in range(3):
    pcart += p[i]*l.matrix[i]/l.lengths[i]
pnorm = np.linalg.norm(pcart)
# print results
print('Polarization (uC/cm2):')
print('crys:',p)
print('quanta:',pq)
print('norm:', pnorm)
