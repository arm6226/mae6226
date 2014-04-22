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
#starting comment goes here
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

pn = np.empty(N,dtype=object)
#eint = np.empty(N,dtype=object)
minf = np.empty((N,N,3,3),dtype=object)

#lid
for k in range(Ni):
    xo,xe=(k)*h,(k+1)*h
    yo,ye=L,L
    pn[k]=Panel(xo,yo,xe,ye,L)
#rhs
for k in range(Ni):
    k1=Ni+k
    xo,xe=L,L
    yo,ye=L-k*h,L-(k+1)*h
    pn[k1]=Panel(xo,yo,xe,ye,L)
#bottom
for k in range(Ni):
    k1=2*Ni+k
    yo,ye=0,0
    xo,xe=L-k*h,L-(k+1)*h
    pn[k1]=Panel(xo,yo,xe,ye,L)
#LHS
for k in range(Ni):
    k1=3*Ni+k
    xo,xe=0,0
    yo,ye=k*h,(k+1)*h
    pn[k1]=Panel(xo,yo,xe,ye,L)
#include the geometry plot here
def I(pj,pk):
    def func(x,y,i,j):
        xh=np.zeros((3,1))
        xh[1]=x-pk.xc
        xh[2]=y-pk.yc
        r=1e-50+sqrt((x-pk.xc)**2+(y-pk.yc)**2)
        fn=0.
        fn=-log(r)*(1-abs(i-j))+xh[i]*xh[j]/r**2
        return fn
    ii=np.zeros((3,3))
    for i in [1,2]:
        for j in [1,2]:
            if (pj.loc=='t'):ii[i,j]=integrate.quad(lambda x:func(x,L,i,j),pj.xa,pj.xb)[0]
            elif (pj.loc=='b'):ii[i,j]=integrate.quad(lambda x:func(x,0.,i,j),pj.xa,pj.xb)[0]
            elif (pj.loc=='l'):ii[i,j]=integrate.quad(lambda y:func(0.,y,i,j),pj.ya,pj.yb)[0]
            elif (pj.loc=='r'):ii[i,j]=integrate.quad(lambda y:func(L,y,i,j),pj.ya,pj.yb)[0]
    return ii
#compute SLP
#for present location, center of ik'th panel
minf=np.zeros((N,N,3,3))
for ik in range(N):
    for i in range(N):
        eint=I(pn[i],pn[ik])
        for m in [1,2]:
            for n in [1,2]:
                minf[ik,k,m,n]=eint[m,n]/(4*pi*mu)

#define DLP integral





#compute DLP
for k in range(Ni):
