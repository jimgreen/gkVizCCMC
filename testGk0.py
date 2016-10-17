import gk0
import matplotlib.pyplot as plt

f = gk0.File('asym_q_7.h5','r')
x = f.getCoordinates('x')
y = f.getCoordinates('y')

rho_e = f.getField(0)
rhovx_e = f.getField(1)
vx_e = rhovx_e / rho_e

plt.pcolormesh(x,y,vx_e)
plt.title('vx_e')

ix=slice(None,512)
iy=slice(256,None)
x_ = f.getCoordinates('x')[ix]
y_ = f.getCoordinates('y')[iy]
Ez_ = f.getField(12, ix=ix,iy=iy)
plt.figure()
plt.pcolormesh(x_,y_,Ez_)
plt.title('Ez in upper left block')

plt.show()
