import h5py
import matplotlib.pyplot as plt

f  = h5py.File('asym_q_7.h5','r')
sg = f['StructGrid']
sgf = f['StructGridField']

rho_e = sgf[:,:,0].T
rhovx_e = sgf[:,:,1].T
vx_e = rhovx_e / rho_e

(xlo,ylo) = sg.attrs['vsLowerBounds']
(xup,yup) = sg.attrs['vsUpperBounds']

plt.imshow(vx_e, origin='lower', extent=[xlo,xup,ylo,yup])
plt.show()
