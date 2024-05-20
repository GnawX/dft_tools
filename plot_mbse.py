import matplotlib
import matplotlib.pyplot as plt
import numpy as np

data1 = np.loadtxt('optics.dat')
data2 = np.loadtxt('mbse/optics_model_exp.dat')
data3 = np.loadtxt('mbse/optics_model_quadratic.dat')
x1 = data1[:,0]
y1 = data1[:,1]
x2 = data2[:,0]
y2 = data2[:,1]
x3 = data3[:,0]
y3 = data3[:,1]

xl = r'$\mathit{h} \nu \mathrm{(eV)}$'
yl = r'$\epsilon_2$'

fig, ax = plt.subplots()
ax.plot(x1,y1,label='RPA')
ax.plot(x2,y2, linestyle='--',label=r'$\mathrm{model}\/\epsilon\/exponential$')
ax.plot(x3,y3, linestyle='-.', label=r'$\mathrm{model}\/\epsilon\/quadratic$')

ax.set(xlabel=xl,ylabel=yl,
       xlim=(5,10),ylim=(0,5))
ax.legend()
plt.show()
