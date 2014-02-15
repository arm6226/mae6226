#Infinite Row of Vortices
import numpy as np
import matplotlib.pyplot as plt
from math import *
plt.close('all')
N = 50                            # Number of points in each direction
xStart,xEnd = -2.0,2.0            # x-direction boundaries
yStart,yEnd = -1.0,1.0            # y-direction boundaries
x = np.linspace(xStart,xEnd,N)    # x 1D-array
y = np.linspace(yStart,yEnd,N)    # y 1D-array
X,Y = np.meshgrid(x,y)            # generation of the mesh grid

gamma = -5.0                  # strength of the vortex
a=0.2
xv=np.linspace(xStart,xEnd,np.abs(xStart-xEnd)/a)
yv=np.zeros_like(xv)

u=gamma/(2*a)*np.sinh(2*pi*Y/a)/(np.cosh(2*pi*Y/a)-np.cos(2*pi*X/a))
v=gamma/(2*a)*np.sin(2*pi*X/a)/(np.cosh(2*pi*Y/a)-np.cos(2*pi*X/a))
#psi=gamma/(2*pi)*np.log(0.5*(np.cosh(2*pi*Y/a)-np.cosh(2*pi*X/a)))

size = 10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.streamplot(X,Y,u,v,\
         density=2,linewidth=1,arrowsize=1,arrowstyle='->')
plt.title('Infinite Row of Vortices', fontsize=16)
plt.scatter(xv,yv,c='#CD2305',s=80,marker='o')
plt.show()