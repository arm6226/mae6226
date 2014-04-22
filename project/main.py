import numpy as np
from scipy import integrate
from math import *
import matplotlib.pyplot as plt
plt.close('all')
Ni = 16;       # Number of elements per side
N  = 4*Ni;    # Total number of elements
L = 1.0;   # square box size
h = L/Ni;  # element length
mu = 1.0;  # Fluid viscosity
U  = 2.0;  # Lid velocity
#declarations

pn = np.empty(N,dtype=object)
minf=np.zeros((N,N,3,3))
dlpsum=np.zeros((3,1))
dd=np.zeros((2*N,1))
A=np.zeros((2*N,2*N))
c1=np.zeros((2*Ni,1))
c2=np.zeros((6*Ni,1))
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
#compute SLP
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
#compute DLP
def DLP(pj,pk):
    dlpint=np.zeros((3,1))
    def func2(x,y,i,j):
        xh=np.zeros((3,1))
        xh[1]=x-pk.xc
        xh[2]=y-pk.yc
        r=1e-50+sqrt((x-pk.xc)**2+(y-pk.yc)**2)
        T = -4*xh[1]*xh[i]*xh[j]/r**4;
        return T
    dlpint[1]=integrate.quad(lambda x:func2(x,L,1,2),pj.xa,pj.xb)[0]
    dlpint[2]=integrate.quad(lambda x:func2(x,L,2,2),pj.xa,pj.xb)[0]
    return dlpint
#for present location, center of ik'th panel
for ik in range(N):
    for k in range(N):
        eint=I(pn[k],pn[ik]) #calculating SLP
        for m in [1,2]:
            for n in [1,2]:
                minf[ik,k,m,n]=eint[m,n]/(4*pi*mu)
#compute DLP
    dlpsum[1]=0.
    dlpsum[2]=0.
    for k in range(Ni):
        dl=DLP(pn[k],pn[ik])# calculate DLP
        dlpsum[1]=dlpsum[1]+(1/(4*pi))*U*dl[1]
        dlpsum[2]=dlpsum[2]+(1/(4*pi))*U*dl[2]
    ik2=2*ik
    dd[ik2]=dlpsum[1]
    dd[ik2+1]=dlpsum[2]

#coeff matrix ----A
for ik in range(N):
    ik2=2*ik
    for k in range(N):
        k2=2*k
        A[ik2,k2]     = minf[ik,k,1,1];
        A[ik2,k2+1]   = minf[ik,k,1,2];
        A[ik2+1,k2]   = minf[ik,k,2,1];
        A[ik2+1,k2+1] = minf[ik,k,2,2];

#RHS---b
for ik in range(Ni):
    ik2=2*ik
    c1[ik2]=1
    c1[ik2+1]=0
for ik in range(3*Ni):
    ik2=2*ik
    c2[ik2]=1
    c2[ik2+1]=0

b= dd - 0.5*U*(np.concatenate([c1,c2]))

#solve the system
var=np.linalg.solve(A,b)

#get the element tractions
fx=var[::2]
fy=var[1::2]

plt.plot(fx,'g-')
plt.plot(fy,'r-')

plt.show()
