from netCDF4 import Dataset
import numpy as np

#model dielectric function
#em1s(q,G,G)=1-(1-em1s000)exp(-|q+G|^2/4 mu^2)
em1s000 = 0.27
mu = 0.61 # screening length bohr^-1
fdb1 = '../SAVE/ns.db1'
fem1sdb = 'ndb.em1s'

# read lattice (cartesian au)
data = Dataset(fdb1)
alat = data['LATTICE_PARAMETER'][:]
data.close()

# read qpoint and gvectors (cartesian)
data = Dataset(fem1sdb)
qpts = data['HEAD_QPT'][:].T
gvecs = data['X_RL_vecs'][:].T
data.close()
nqpts = len(qpts)
ngvecs = len(gvecs)


# write v\chi(G,G',q)
for nq in range(nqpts):
    print('creating fragment_%d'%(nq+1))
    filename = "%s_fragment_%d"%(fem1sdb,nq+1)
    data = Dataset(filename,'w')
    dimg = 'D_{:0=10d}'.format(ngvecs)
    dim1 = 'D_{:0=10d}'.format(1)
    dim2 = 'D_{:0=10d}'.format(2)
    dim6 = 'D_{:0=10d}'.format(6)
    data.createDimension(dim1, 1)
    data.createDimension(dim2, 2)
    data.createDimension(dim6, 6)
    data.createDimension(dimg, ngvecs)
    pfq = data.createVariable('FREQ_PARS_sec_iq%d'%(nq+1),'f4',(dim6,))
    fq = data.createVariable('FREQ_sec_iq%d'%(nq+1),'f4',(dim1,dim6,))
    eps = data.createVariable('X_Q_%d'%(nq+1),'f4',(dim1,dimg,dimg,dim2,))
    pfq[:] = 0.
    pfq[0] = 1.
    pfq[3] = 3.6749327e-05
    pfq[4] = 3.6749327e-05
    pfq[5] = 1.0
    fq[:,:] = 0.
    fq[0,1] = 3.6749327e-05
    eps[:,:,:,:] = 0.
    for ng in range(ngvecs):
        qplusg = np.linalg.norm((qpts[nq]+gvecs[ng])/alat*2.0*np.pi)
        eps[0,ng,ng,0] = (em1s000-1.0)*np.exp(-qplusg**2/4.0/mu**2) 
    data.close()
