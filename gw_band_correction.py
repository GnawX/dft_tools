import numpy as np
from pymatgen.io.vasp.outputs import Eigenval
from pymatgen.electronic_structure.core import Spin
import matplotlib
import matplotlib.pyplot as plt
import scipy.optimize as opt


vbcut = -4
cbcut = 4

gw = Eigenval('EIGENVAL_GW')
dft = Eigenval('EIGENVAL_DFT')

nk = dft.nkpt
nb = dft.nbands
ne = dft.nelect
ivbm = int(np.rint(ne/2))

vbm1 = dft.eigenvalue_band_properties[2]
vbm2 = gw.eigenvalue_band_properties[2]
cbm1 = dft.eigenvalue_band_properties[1]
cbm2 = gw.eigenvalue_band_properties[1]

vb_dft = dft.eigenvalues[Spin.up][:,:ivbm,0].flatten()-vbm1
vb_gw = gw.eigenvalues[Spin.up][:,:ivbm,0].flatten()-vbm2
cb_dft = dft.eigenvalues[Spin.up][:,ivbm:,0].flatten()-cbm1
cb_gw = gw.eigenvalues[Spin.up][:,ivbm:,0].flatten()-cbm2

vb = np.column_stack((vb_dft, vb_gw))
cb = np.column_stack((cb_dft, cb_gw))

#np.savetxt('valence_dft_gw.dat',vb,fmt='%10.4f')
#np.savetxt('conduction_dft_gw.dat',cb,fmt='%10.4f')

maskv = vb_dft > vbcut
maskc = cb_dft < cbcut

x1 = vb_dft[maskv]
y1 = vb_gw[maskv]
x2 = cb_dft[maskc]
y2 = cb_gw[maskc]

xl = r'DFT (eV)'
yl = r'GW (eV)'
xv = np.linspace(vbcut,0,1000)
xc = np.linspace(0,cbcut,1000)

def func(x,a):
    return a*x

popt1, pcov1 = opt.curve_fit(func, x1, y1)
popt2, pcov2 = opt.curve_fit(func, x2, y2)
print(popt1, popt2)

plt.scatter(x1,y1,c='C0',marker='+',label='VB')
plt.scatter(x2,y2,c='C1',marker='x',label='CB')
plt.plot(xv,func(xv,*popt1),'b-', label='fit-VB')
plt.plot(xc,func(xc,*popt2),'r-', label='fit-CB')

plt.xlabel(xl)
plt.ylabel(yl)
plt.legend()
plt.show()
