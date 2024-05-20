import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt

data1 = np.loadtxt('eps_G.dat',usecols=(0,1))
eps = 1./data1[0,1]
x1 = data1[:,0]
y1 = data1[:,1]
xm = max(x1)
xl = r'$|q+G|\/(\AA^{-1})$'
yl = r'$\epsilon^{-1}$'
xf = np.linspace(0.,xm,1000)

def func1(x,a):
    return 1.-(1.-1./eps)*np.exp(-x**2/4./a**2)
def func2(x,a,b):
    return 1./(1. + 1./(1./(eps-1) + a*x**2 + b*x**4))

popt1, pcov1 = opt.curve_fit(func1, x1, y1)
popt2, pcov2 = opt.curve_fit(func2, x1, y1)
print(popt1)
print(popt2[0]*4*np.pi**2,popt2[1]*(2*np.pi)**4)

#fig, ax = plt.subplots()
#ax.plot(x1,y1,'ko',markersize=5, label='RPA')
plt.scatter(x1,y1,c='k',label='RPA')
plt.plot(xf,func1(xf,*popt1),'b-', label='fit-exponential')
plt.plot(xf,func2(xf,*popt2),'r-', label='fit-quadratic')

plt.xlabel(xl)
plt.ylabel(yl)
       #xlim=(5,10),ylim=(0,5))
plt.legend()
plt.show()
