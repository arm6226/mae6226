import numpy as np
import matplotlib.pyplot as plt
from math import *
plt.close('all')
N = 100                           # Number of points in each direction
xStart,xEnd = -1.0,1.0            # x-direction boundaries
yStart,yEnd = -1.5,1.5            # y-direction boundaries
x = np.linspace(xStart,xEnd,N)    # x 1D-array
y = np.linspace(yStart,yEnd,N)    # y 1D-array
X,Y = np.meshgrid(x,y)            # generation of the mesh grid

Uinf = 1.0

uFreestream = Uinf*np.ones((N,N),dtype=float)
vFreestream = np.zeros((N,N),dtype=float)

class Source:
    def __init__(self,strength,x,y):
        self.strength = strength
        self.x,self.y = x,y
    # get the velocity field
    def velocity(self,X,Y):
        self.u = self.strength/(2*pi)*(X-self.x)/((X-self.x)**2+(Y-self.y)**2)
        self.v = self.strength/(2*pi)*(Y-self.y)/((X-self.x)**2+(Y-self.y)**2)
    # get the stream-function
    def streamFunction(self,X,Y):
        self.psi = self.strength/(2*pi)*np.arctan2((Y-self.y),(X-self.x))
        
# definition of my sources
Nsources = 11
strength = 5.0
strengthSource = strength/Nsources
xSource = 0.0
ySource = np.linspace(-1.0,1.0,Nsources)

# creation of the source line
sources = np.empty( Nsources, dtype=object)

for i in range(Nsources):
    sources[i] = Source(strengthSource,xSource,ySource[i])
    sources[i].velocity(X,Y)

# superposition
u = uFreestream.copy()
v = vFreestream.copy()

for s in sources:
    u = np.add(u,s.u)
    v = np.add(v,s.v)
    
# plotting

size = 6
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.grid(True)
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)

plt.streamplot( X,Y, u,v, density=2,linewidth=1,arrowsize=1,arrowstyle='->')
plt.scatter( xSource*np.ones(Nsources,dtype=float), ySource,c='#CD2305',s=80,marker='o')
velocity = plt.contourf( X,Y,np.sqrt(u**2+v**2),levels=np.linspace(0.0,0.1,10))
cbar = plt.colorbar(velocity, ticks=[0, 0.05, 0.1], orientation='horizontal')
cbar.set_label('Velocity magnitude',fontsize=16);

from scipy import integrate

sigma = 0.5    # strength of the source sheet
uPanel    = np.empty((N,N),dtype=float)
vPanel    = np.empty((N,N),dtype=float)

# boundaries of the sheet
ymin,ymax = -1.0,1.0

# computing the velocity field
for i in range(N):
    for j in range(N):
        
        integrand = lambda s : X[i,j]/(X[i,j]**2+(Y[i,j]-s)**2)
        uPanel[i,j] = sigma/(2*pi)*integrate.quad( integrand, ymin, ymax)[0]
        
        integrand = lambda s: (Y[i,j]-s)/(X[i,j]**2+(Y[i,j]-s)**2)
        vPanel[i,j] = sigma/(2*pi)*integrate.quad( integrand, ymin, ymax)[0]

# superposition to the uniform flow
u2   = np.add(uFreestream,uPanel)
v2   = np.add(vFreestream,vPanel)

# plotting


size = 8
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.grid(True)
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.streamplot( X,Y, u2,v2, density=2,linewidth=1,arrowsize=1,arrowstyle='->')
plt.axvline( 0.0,(ymin-yStart)/(yEnd-yStart),(ymax-yStart)/(yEnd-yStart),\
            color='#CD2305',linewidth=4)
velocity = plt.contourf( X,Y, np.sqrt(u2**2+v2**2),levels=np.linspace(0.0,0.1,10))
cbar = plt.colorbar(velocity, ticks=[0, 0.05, 0.1], orientation='horizontal')
cbar.set_label('Velocity magnitude',fontsize=16);
    
plt.show()