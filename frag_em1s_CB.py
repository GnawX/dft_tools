from netCDF4 import Dataset
import numpy as np

#model dielectric function
#em1s(q,G,G)=1/(1+1/(1/(1/em1s000-1) + a*|q+G|^2 + b|q+G|^4))
em1s000 = 2.624395489692687988e-01
a,b = 1.5092446601249485,1.1383845995046515
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
        eps[0,ng,ng,0] = 1./(1. + 1./(1./(1./em1s000-1) + a*qplusg**2 + b*qplusg**4)) - 1.0
    data.close()
