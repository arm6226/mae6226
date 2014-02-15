#Finite Row of Vortices
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
xVortex,yVortex = 0.0,0.0    # location of the vortex

def getStreamFunctionVortex(strength,xv,yv,X,Y):
    psi = strength/(4*pi)*np.log((X-xv)**2+(Y-yv)**2)
    return psi

def getVelocityVortex(strength,xv,yv,X,Y):
    u = + strength/(2*pi)*(Y-yv)/((X-xv)**2+(Y-yv)**2)
    v = - strength/(2*pi)*(X-xv)/((X-xv)**2+(Y-yv)**2)
    return u,v

uVortex=np.zeros_like(X)
vVortex=np.zeros_like(X)
psiVortex=np.zeros_like(X)
        
#number of point-vortices
K=10
#x coordinate of the vortices
xv=np.linspace(xStart,xEnd,K)
yv=0

# computing the velocity components on the mesh grid
for i, xVort in enumerate(xv):
    uVortex1,vVortex1 = getVelocityVortex(gamma,xVort,yv,X,Y)
    uVortex=uVortex+uVortex1
    vVortex=vVortex+vVortex1
# computing the stream-function on the mesh grid
    psiVortex = psiVortex+getStreamFunctionVortex(gamma,xVort,yv,X,Y)
#    psiVortex=psiVortex+psiVortex1

yvp=np.zeros_like(xv)
# plotting
size = 10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.streamplot(X,Y,uVortex,vVortex,\
               density=3.5,linewidth=1,arrowsize=1,arrowstyle='->')
plt.scatter(xv,yvp,c='#CD2305',s=80,marker='o')
plt.title('Finite Row of Vortices', fontsize=16)
plt.show()