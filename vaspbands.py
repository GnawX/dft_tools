import py4vasp

calc = py4vasp.Calculation.from_path('.')
bands = calc.band.to_dict(source='kpoints_opt')['bands']
kpts = calc.band.to_dict(source='kpoints_opt')['kpoint_distances']
nk, nb = bands.shape

with open('bands.dat','w') as f:
     for ib in range(nb):
         f.write('\n')
         for ik in range(nk):
             f.write(('{:12.4f}'+'{:12.4f}\n').format(kpts[ik],bands[ik,ib])) 
