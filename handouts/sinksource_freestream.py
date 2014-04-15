import numpy as np
import matplotlib.pyplot as plt
from math import *

plt.close('all')

N = 200                           # Number of points in each direction
xStart,xEnd = -4.0,4.0            # x-direction boundaries
yStart,yEnd = -2.0,2.0            # y-direction boundaries
x = np.linspace(xStart,xEnd,N)    # x 1D-array
y = np.linspace(yStart,yEnd,N)    # y 1D-array
X,Y = np.meshgrid(x,y)            # generation of the mesh grid

Uinf = 1.0        # freestream speed

# computing the velocity components on the mesh grid
uFreestream = Uinf*np.ones((N,N),dtype=float)
vFreestream = np.zeros((N,N),dtype=float)

# computing the stream-function on the mesh grid
psiFreestream = Uinf*Y

# function to compute the velocity field of a source/sink
def getVelocity(strength,xs,ys,X,Y):
    u = strength/(2*pi)*(X-xs)/((X-xs)**2+(Y-ys)**2)
    v = strength/(2*pi)*(Y-ys)/((X-xs)**2+(Y-ys)**2)
    return u,v

# function to compute the stream-function of a source/sink
def getStreamFunction(strength,xs,ys,X,Y):
    psi = strength/(2*pi)*np.arctan2((Y-ys),(X-xs))
    return psi

strengthSource = 5.0         # strength of the source
xSource,ySource = -1.0,0.0   # location of the source

# computing the velocity components
uSource,vSource = getVelocity(strengthSource,xSource,ySource,X,Y)

# computing the stream-function
psiSource = getStreamFunction(strengthSource,xSource,ySource,X,Y)

# superposition of the source on the freestream
u = uFreestream + uSource
v = vFreestream + vSource
psi = psiFreestream + psiSource

size = 10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.grid(True)
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.streamplot(X,Y,u,v,density=2.0,linewidth=1,arrowsize=1,arrowstyle='->')
plt.scatter(xSource,ySource,c='#CD2305',s=80,marker='o')

# computing the stagnation point
xStagnation = xSource - strengthSource/(2*pi*Uinf)
yStagnation = ySource

# adding the stagnation point to the figure
plt.scatter(xStagnation,yStagnation,c='g',s=80,marker='o')

# adding the dividing line to the figure
plt.contour(X,Y,psi,\
            levels=[-strengthSource/2,+strengthSource/2],\
            colors='#CD2305',linewidths=2,linestyles='solid');


strengthSink = -5.0     # strength of the sink
xSink,ySink = 1.0,0.0   # location of the sink

# computing the velocity field on the mesh grid
uSink,vSink = getVelocity(strengthSink,xSink,ySink,X,Y)

# computing the stream-function on the grid mesh
psiSink = getStreamFunction(strengthSink,xSink,ySink,X,Y)


# superposition of a source and a sink on the freestream
u = uFreestream + uSource + uSink
v = vFreestream + vSource + vSink
psi = psiFreestream + psiSource + psiSink

# plotting
size = 10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.streamplot(X,Y,u,v,density=2.0,linewidth=1,arrowsize=1,arrowstyle='->')
plt.scatter([xSource,xSink],[ySource,ySink],c='#CD2305',s=80,marker='o')
plt.contour(X,Y,psi,levels=[0.0],colors='#CD2305',linewidths=2,linestyles='solid');

 # computing the pressure coefficient
Cp = 1.0-(u**2+v**2)/Uinf**2

# plotting
size = 10
plt.figure(figsize=(1.1*size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
contf = plt.contourf(X,Y,Cp,levels=np.linspace(-2.0,1.0,100),extend='both')
cbar = plt.colorbar(contf)
cbar.set_label(r'$C_p$',fontsize=16)
cbar.set_ticks([-2.0,-1.0,0.0,1.0])
plt.scatter([xSource,xSink],[ySource,ySink],c='#CD2305',s=80,marker='o')
plt.contour(X,Y,psi,\
            levels=[0.0],\
            colors='#CD2305',linewidths=2,linestyles='solid');
plt.show()
            
    