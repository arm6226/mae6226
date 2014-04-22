import numpy as np
from scipy import integrate
from math import *
import matplotlib.pyplot as plt
plt.close('all')
Ni = 8;       # Number of elements per side
N  = 4*Ni;    # Total number of elements

L = 1.0;   # square box size

h = L/Ni;  # element length

mu = 1.0;  # Fluid viscosity
U  = 2.0;  # Lid velocity

####
# Number elements starting from
# left end of lid, and assign
# coordinates
#
# elmx(k,1/2) : x coord of start (0) or end (1)
#                 of element k
#
# elmy(k,1/2) : y coord of start (0) or end (1)
#                 of element k
####
#                U
# ______________--->_______________
#|              1                 |
#|                                |
#|                                |
#|                                |
#|                             2  |
#|4                               |
#|                                |
#|                                |
#|               3                |
#|________________________________|
#                 L
class Panel:
    def __init__(self,xa,ya,xb,yb,L):
        self.L=L
        self.xa,self.ya = xa,ya                     # 1st end-point
        self.xb,self.yb = xb,yb                     # 2nd end-point
        self.xc,self.yc = (xa+xb)/2,(ya+yb)/2       # control point
        if(abs(self.xc)<=1.e-5): self.loc='l'
        if(abs(self.yc)<=1.e-5): self.loc='b'
        if(abs(self.xc-self.L)<=1.e-5): self.loc='r'
        if(abs(self.yc-self.L)<=1.e-5): self.loc='t'
        if(self.loc=='t' or self.loc=='b'):
            self.ielm=0
        else:self.ielm=1

panel = np.empty(N,dtype=object)
#lid
for k in range(Ni):
    xo=(k)*h
    xe=(k+1)*h
    yo=L
    ye=L
    panel[k]=Panel(xo,yo,xe,ye,L)

#rhs
for k in range(Ni):
    k1=Ni+k
    xo=L
    xe=L
    yo=L-k*h
    ye=L-(k+1)*h;
    panel[k1]=Panel(xo,yo,xe,ye,L)
#bottom
for k in range(Ni):
    k1=2*Ni+k
    yo=0
    ye=0
    xo=L-k*h
    xe=L-(k+1)*h
    panel[k1]=Panel(xo,yo,xe,ye,L)
#LHS
for k in range(Ni):
    k1=3*Ni+k
    xo=0
    xe=0
    yo=k*h
    ye=(k+1)*h
    panel[k1]=Panel(xo,yo,xe,ye,L)


valX,valY = 0.2,0.2
xmin,xmax = min([p.xa for p in panel]),max([p.xa for p in panel])
ymin,ymax = min([p.ya for p in panel]),max([p.ya for p in panel])
xStart,xEnd = xmin-valX*(xmax-xmin),xmax+valX*(xmax-xmin)
yStart,yEnd = ymin-valY*(ymax-ymin),ymax+valY*(ymax-ymin)
size = 10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))

plt.grid(True)
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
#plt.plot(xp,yp,'k-',linewidth=2)
plt.plot(np.append([p.xa for p in panel],panel[0].xa),\
        np.append([p.ya for p in panel],panel[0].ya),\
        linestyle='-',linewidth=1,\
        marker='o',markersize=6,color='#CD2305');

plt.plot([p.xc for p in panel if p.loc=='l'],\
		[p.yc for p in panel if p.loc=='l'],\
		'ko',linewidth=5)
plt.plot([p.xc for p in panel if p.loc=='r'],\
		[p.yc for p in panel if p.loc=='r'],\
		'ko',linewidth=5)
plt.plot([p.xc for p in panel if p.loc=='t'],\
		[p.yc for p in panel if p.loc=='t'],\
		'co',linewidth=5)
plt.plot([p.xc for p in panel if p.loc=='b'],\
		[p.yc for p in panel if p.loc=='b'],\
		'ko',linewidth=5)


plt.show()
