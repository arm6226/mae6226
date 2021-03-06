import numpy as np
from scipy import integrate
from math import *
import matplotlib.pyplot as plt
plt.close('all')
Ni  =  32;       # Number of elements per side
N   =  4*Ni;    # Total number of elements
Lx  =  1.0;   # square box size
Ly  =  1.0;
hx  =  Lx/Ni;  # element length
hy  =  Ly/Ni;
mu  =  1.0;  # Fluid viscosity
U   =  2.0;  # Lid velocity
#declarations

pn = np.empty(N,dtype=object)
minf=np.zeros((N,N,3,3))
dlpsum=np.zeros((3,1))
dd=np.zeros((2*N,1))
A=np.zeros((2*N,2*N))
c1=np.zeros((2*Ni,1))
c2=np.zeros((6*Ni,1))

#define the class panel
class Panel:
    def __init__(self,xa,ya,xb,yb,Lx,Ly):
        self.Lx,self.Ly = Lx,Ly
        self.xa,self.ya = xa,ya                     # 1st end-point
        self.xb,self.yb = xb,yb                     # 2nd end-point
        self.xc,self.yc = (xa+xb)/2,(ya+yb)/2       # control point
        # orientation of the panel
        #location of the element,
        if(abs(self.xc)<=1.e-8): self.loc='l'
        if(abs(self.yc)<=1.e-8): self.loc='b'
        if(abs(self.xc-self.Lx)<=1.e-8): self.loc='r'
        if(abs(self.yc-self.Ly)<=1.e-8): self.loc='t'
        #cartesian components of traction vector
        self.fx=0.
        self.fy=0.
#geometry
for i in range(N):
    if(i<Ni):
        pn[i]=Panel(i*hx,Ly,(i+1)*hx,Ly,Lx,Ly)
    elif(Ni<=i<2*Ni):
        ii=i-Ni
        pn[i]=Panel(Lx,(Ly-ii*hy),Lx,(Ly-(ii+1)*hy),Lx,Ly)
    elif (2*Ni<=i<3*Ni):
        ii=i-2*Ni
        pn[i]=Panel(Lx-ii*hx,0.,Lx-(ii+1)*hx,0.,Lx,Ly)
    elif (3*Ni<=i<4*Ni):
        ii=i-3*Ni
        pn[i]=Panel(0.,(ii*hy),0.,((ii+1)*hy),Lx,Ly)

##include the geometry plot here
#valX,valY = 0.2,0.2
#xmin,xmax = min([p.xa for p in pn]),max([p.xa for p in pn])
#ymin,ymax = min([p.ya for p in pn]),max([p.ya for p in pn])
#xStart,xEnd = xmin-valX*(xmax-xmin),xmax+valX*(xmax-xmin)
#yStart,yEnd = ymin-valY*(ymax-ymin),ymax+valY*(ymax-ymin)
#size = 10
#plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
#
#plt.grid(True)
#plt.xlabel('x',fontsize=16)
#plt.ylabel('y',fontsize=16)
#plt.xlim(xStart,xEnd)
#plt.ylim(yStart,yEnd)
#
#plt.plot(np.append([p.xa for p in pn],pn[0].xa),\
#        np.append([p.ya for p in pn],pn[0].ya),\
#        linestyle='-',linewidth=1,\
#        marker='o',markersize=6,color='#CD2305');
#
#plt.plot([p.xc for p in pn],[p.yc for p in pn],'ko',linewidth=5)
#
#plt.axis('equal')
#compute SLP
def SLP(pj,x0,y0):
    def func(x,y,i,j):
        xh=np.zeros((3,1))
        xh[1]=x-x0
        xh[2]=y-y0
        r=1e-50+np.sqrt((x-x0)**2+(y-y0)**2)
        fn=0.
        fn=-log(r)*(1-abs(i-j))+xh[i]*xh[j]/r**2
        return fn
    ii=np.zeros((3,3))
    for i in [1,2]:
        for j in [1,2]:
            if (pj.loc=='t'):ii[i,j]=integrate.quad(lambda x:func(x,Ly,i,j),pj.xa,pj.xb)[0]
            elif (pj.loc=='b'):ii[i,j]=integrate.quad(lambda x:func(x,0.,i,j),pj.xa,pj.xb)[0]
            elif (pj.loc=='l'):ii[i,j]=integrate.quad(lambda y:func(0.,y,i,j),pj.ya,pj.yb)[0]
            elif (pj.loc=='r'):ii[i,j]=integrate.quad(lambda y:func(Lx,y,i,j),pj.ya,pj.yb)[0]
    return ii
#compute DLP
def DLP(pj,x0,y0):
    dlpint=np.zeros((3,1))
    def func2(x,y,i,j):
        xh=np.zeros((3,1))
        xh[1]=x-x0
        xh[2]=y-y0
        r=1e-50+sqrt((x-x0)**2+(y-y0)**2)
        T = -4*xh[1]*xh[i]*xh[j]/r**4;
        return -T #normal pointing in -y dir
    dlpint[1]=integrate.quad(lambda x:func2(x,Ly,1,2),pj.xa,pj.xb)[0]
    dlpint[2]=integrate.quad(lambda x:func2(x,Ly,2,2),pj.xa,pj.xb)[0]
    return dlpint
#for present location, center of ik'th panel
for ik in range(N):
    for k in range(N):
        eint=SLP(pn[k],pn[ik].xc,pn[ik].yc) #calculating SLP
        for m in [1,2]:
            for n in [1,2]:
                minf[ik,k,m,n]=eint[m,n]/(4*pi*mu)
#compute DLP
    dlpsum[1]=0.
    dlpsum[2]=0.
    for k in range(Ni):
        dl=DLP(pn[k],pn[ik].xc,pn[ik].yc) # calculate DLP
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
    c2[ik2]=0
    c2[ik2+1]=0

b= dd - 0.5*U*(np.concatenate([c1,c2]))

#solve the system
var=np.linalg.solve(A,b)

#get the element tractions and put them in the panels
fx=var[::2]
fy=var[1::2]
for i in range(N):
    pn[i].fx=fx[i]
    pn[i].fy=fy[i]

#function to calculate velocity over a meshgrid
def VelocityField(p,X,Y):
    Nx,Ny=X.shape
    N=p.size
    u,v = np.empty((Nx,Ny),dtype=float),np.empty((Nx,Ny),dtype=float)
    for i in range(Nx):
        for j in range(Ny):
            ux=0.
            uy=0.
            for k in range(N):
                eint=SLP(p[k],X[i,j],Y[i,j])
                ux = ux - (1/(4*pi*mu))*(eint[1,1]*p[k].fx + eint[2,1]*p[k].fy);
                uy = uy - (1/(4*pi*mu))*(eint[1,2]*p[k].fx + eint[2,2]*p[k].fy);
            dlpsum[1]=0.
            dlpsum[2]=0.
            for k in range(N/4):
                dl=DLP(p[k],X[i,j],Y[i,j])
                dlpsum[1]=dlpsum[1]+(1/(4*pi))*U*dl[1]
                dlpsum[2]=dlpsum[2]+(1/(4*pi))*U*dl[2]
            ux = ux + dlpsum[1]
            uy = uy + dlpsum[2]
            u[i,j]=ux
            v[i,j]=uy
    return u,v

#definition of meshgrid
Nx,Ny=20,20
X,Y=np.meshgrid(np.linspace(0.,0.01*Lx,Nx),np.linspace(0.,0.01*Ly,Ny))

u,v=VelocityField(pn,X,Y)

valX,valY = 0.1,0.1
xmin,xmax = min([p.xa for p in pn]),max([p.xa for p in pn])
ymin,ymax = min([p.ya for p in pn]),max([p.ya for p in pn])
xStart,xEnd = xmin-valX*(xmax-xmin),xmax+valX*(xmax-xmin)
yStart,yEnd = ymin-valY*(ymax-ymin),ymax+valY*(ymax-ymin)
size=12
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.streamplot(X,Y,u,v,density=10,linewidth=1,arrowsize=1,arrowstyle='->')
#plt.quiver(X,Y,u,v)
#plt.plot(np.append([p.xa for p in pn],pn[0].xa),\
#        np.append([p.ya for p in pn],pn[0].ya),\
#        linestyle='-',linewidth=1,\
#        marker='o',markersize=6,color='#CD2305');
plt.xlim(0.,0.011)
plt.ylim(0.,0.011)
plt.axis('equal')
plt.title('Contour of velocity field');
plt.show()
