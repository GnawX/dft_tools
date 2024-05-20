from netCDF4 import Dataset
import numpy as np

job = './'
fdb1 = '../SAVE/ns.db1'
fem1sdb = job+'ndb.em1s'

# read lattice (cartesian au)
data = Dataset(fdb1)
a = data['LATTICE_VECTORS'][:].T
alat = data['LATTICE_PARAMETER'][:]
data.close()
volume = np.linalg.det(a)
#b = np.linalg.inv(a.T)*np.pi*2.0

# read qpoint and gvectors (cartesian)
data = Dataset(fem1sdb)
qpts = data['HEAD_QPT'][:].T
gvecs = data['X_RL_vecs'][:].T
data.close()
nqpts = len(qpts)
ngvecs = len(gvecs)

em1s = np.zeros([nqpts,ngvecs,ngvecs],dtype=np.complex64)
xq = np.zeros([nqpts])
# read v\chi(G,G',q)
for nq in range(nqpts):
    filename = "%s_fragment_%d"%(fem1sdb,nq+1)
    data = Dataset(filename)
    re = data['X_Q_%d'%(nq+1)][0,:,:,0]
    im = data['X_Q_%d'%(nq+1)][0,:,:,1]
    em1s[nq,:,:] = re + 1j*im
    data.close()
    #xq[nq] = np.linalg.norm(qpts[nq]/alat*2.0*np.pi)
em1s += 1.0

#head
#head = np.concatenate((xq,em1s[:,0,0].real),axis=0)
#np.savetxt('em1s_00.dat',np.reshape(head,(nqpts,2),order='F'))

#diagonal
diag = np.zeros([ngvecs*nqpts,2])
count = 0
for nq in range(nqpts):
    for ng in range(ngvecs):
        diag[count,0] = np.linalg.norm((qpts[nq]+gvecs[ng])/alat*2.0*np.pi)
        diag[count,1] = em1s[nq,ng,ng].real
        count += 1 
np.savetxt('em1s_GG.dat',diag)
#
#print('em1s000:',em1s[0,0,0].real)
