import numpy as np
from pymatgen.io.vasp.outputs import Outcar
from pymatgen import Structure
from pymatgen.analysis.ferroelectricity.polarization import Polarization, EnergyTrend

nimages = 21
outcar = []
atoms = []
energy = []

for i in range(nimages):
    fname1 = 'image'+str(i)+'/POSCAR'
    fname2 = 'image'+str(i)+'/OUTCAR'
    atoms += [ Structure.from_file(fname1) ]
    o = Outcar(fname2)
    energy += [o.final_energy]
    outcar += [ o ]
pol = Polarization.from_outcars_and_structures(outcar,atoms,calc_ionic_from_zval=True)

p = pol.get_polarization_change_norm()
pb = pol.get_same_branch_polarization_data()
pq = pol.get_lattice_quanta()
#spa,spb,spc=pol.same_branch_splines()

esp = EnergyTrend(energy).spline()

np.savetxt('pbrach.dat',pb)
np.savetxt('pquanta.dat',pq)
energy = np.asarray(energy)
np.savetxt('energy.dat',np.transpose(energy))
#sp = np.zeros([21,3])
e = np.zeros(21)
for i in range(21):
    e[i] = esp(i)
#    sp[i,0]= spa(i)
#    sp[i,1]= spb(i)
#    sp[i,2]= spc(i)
#np.savetxt('pspline.dat',sp)
np.savetxt('espline.dat',e)
print(p)
