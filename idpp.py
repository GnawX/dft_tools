from ase import io
from ase.neb import NEB

ni = 11
inital = io.read('POSCAR_ferro')
final = io.read('POSCAR_para')

images = [inital]

for i in range(ni):
    images.append(inital.copy())
images.append(final)

neb = NEB(images)
neb.interpolate('idpp')

for i in range(1,ni+1):
    fil = 'POSCAR_' + str(i)
    io.write(fil,images[i],direct=True)
