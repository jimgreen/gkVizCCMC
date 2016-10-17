import gk3
import matplotlib.pyplot as plt

f = gk3.File('asym_q_7.h5')
x = f.getCoordinates('x')
y = f.getCoordinates('y')

vx_e = f.getField('ux_e')
plt.pcolormesh(x,y,vx_e)
plt.title('vx_e as a pre-defined quantity')

vx_e = f.getField('rhoux_e/rho_e')
plt.figure()
plt.pcolormesh(x,y,vx_e)
plt.title('vx_e using calculator on the fly')

ix=slice(None,512)
iy=slice(256,None)
x_ = f.getCoordinates('x')[ix]
y_ = f.getCoordinates('y')[iy]
Ez_ = f.getField('Ez', ix=ix,iy=iy)
plt.figure()
plt.pcolormesh(x_,y_,Ez_)
plt.title('Ez in left upper block')

plt.show()
