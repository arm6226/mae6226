import numpy as np
from scipy import integrate
from math import *
import matplotlib.pyplot as plt
plt.close('all')
Ni  =  32;       # Number of elements per side
N   =  5*Ni;    # Total number of elements
Lx,Ly  =  [1.0,1.0];   # square box size
hx,hy =  [Lx/Ni,Ly/Ni];  # element length
mu  =  1.0;  # Fluid viscosity
U   =  2.0;  # Lid velocity
pn = np.empty(N,dtype=object)
minf=np.zeros((N,N,3,3))
dlpsum=np.zeros((3,1))
dd=np.zeros((2*N,1))
A=np.zeros((2*N,2*N))
c1=np.zeros((2*Ni,1))
c2=np.zeros((8*Ni,1))
class Panel:
    def __init__(self,xa,ya,xb,yb):
        self.xa,self.ya = xa,ya                     # 1st end-point
        self.xb,self.yb = xb,yb                     # 2nd end-point
        self.xc,self.yc = (xa+xb)/2,(ya+yb)/2       # control point
        self.length = sqrt((xb-xa)**2+(yb-ya)**2)
        # orientation of the panel
        if (xb-xa>=0.):
            self.beta=acos((ya-yb)/self.length)
        elif (xb-xa<0.):
            self.beta = pi+acos((ya-yb)/self.length)
        #cartesian components of traction vector
        self.fx=0.
        self.fy=0.
#define the geometry
for i in range(N):
    if(i<Ni):
        pn[i]=Panel(i*hx,Ly,(i+1)*hx,Ly)
    elif(Ni<=i<2*Ni):
        ii=i-Ni
        pn[i]=Panel(Lx,(Ly-ii*hy),Lx,(Ly-(ii+1)*hy))
    elif (2*Ni<=i<3*Ni):
        ii=i-2*Ni
        pn[i]=Panel(Lx-ii*hx,0.,Lx-(ii+1)*hx,0.)
    elif (3*Ni<=i<4*Ni):
        ii=i-3*Ni
        pn[i]=Panel(0.,(ii*hy),0.,((ii+1)*hy))
    elif (4*Ni<=i<5*Ni):
        ii=i-3*Ni
        xa=0.5 +0.2*cos(2*pi*float(Ni-ii)/Ni)
        ya=0.5 +0.2*sin(2*pi*float(Ni-ii)/Ni)
        xb=0.5 +0.2*cos(2*pi*float(Ni-ii-1)/Ni)
        yb=0.5 +0.2*sin(2*pi*float(Ni-ii-1)/Ni)
        pn[i]=Panel(xa,ya,xb,yb)
#display the geometry
valX,valY = 0.1*Lx,0.1*Ly
xmin,xmax = min([p.xa for p in pn]),max([p.xa for p in pn])
ymin,ymax = min([p.ya for p in pn]),max([p.ya for p in pn])
xStart,xEnd = xmin-valX*(xmax-xmin),xmax+valX*(xmax-xmin)
yStart,yEnd = ymin-valY*(ymax-ymin),ymax+valY*(ymax-ymin)
size = 10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.grid(True)
plt.xlabel('x',fontsize=16);plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd);plt.ylim(yStart,yEnd)
plt.plot(np.append([p.xa for p in pn],pn[0].xa),\
        np.append([p.ya for p in pn],pn[0].ya),\
        'ko-');
plt.plot([p.xc for p in pn ],[p.yc for p in pn],'co',linewidth=5)
#compute SLP
def SLP(pj,x0,y0):
    def func(s,i,j):
        xh=np.zeros((3,1))
        x=pj.xa+sin(pj.beta)*s
        y=pj.ya-cos(pj.beta)*s
        xh[1]=x-x0
        xh[2]=y-y0
        r=1e-50+np.sqrt((x-x0)**2+(y-y0)**2)
        fn=-log(r)*(1-abs(i-j))+xh[i]*xh[j]/r**2
        return fn
    ii=np.zeros((3,3))
    for i in [1,2]:
        for j in [1,2]:
            ii[i,j]=integrate.quad(lambda s:func(s,i,j),0.,pj.length)[0]
    return ii
def DLP(pj,x0,y0):
    dlpint=np.zeros((3,1))
    def func2(s,i,j):
        xh=np.zeros((3,1))
        x=pj.xa+sin(pj.beta)*s
        y=pj.ya-cos(pj.beta)*s
        xh[1]=x-x0
        xh[2]=y-y0
        r=1e-50+sqrt((x-x0)**2+(y-y0)**2)
        T = -4*xh[1]*xh[i]*xh[j]/r**4;
        return -T #normal pointing in -y dir
    dlpint[1]=integrate.quad(lambda s:func2(s,1,2),0.,pj.length)[0]
    dlpint[2]=integrate.quad(lambda s:func2(s,2,2),0.,pj.length)[0]
    return dlpint
for ik in range(N):
    for k in range(N):
        eint=SLP(pn[k],pn[ik].xc,pn[ik].yc) #calculating SLP
        for m in [1,2]:
            for n in [1,2]:
                minf[ik,k,m,n]=eint[m,n]/(4*pi*mu)

    dlpsum[1]=0.
    dlpsum[2]=0.
    for k in range(Ni):
        dl=DLP(pn[k],pn[ik].xc,pn[ik].yc) # calculate DLP
        dlpsum[1]=dlpsum[1]+(1/(4*pi))*U*dl[1]
        dlpsum[2]=dlpsum[2]+(1/(4*pi))*U*dl[2]
    ik2=2*ik
    dd[ik2]=dlpsum[1]
    dd[ik2+1]=dlpsum[2]

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
Nx,Ny=32,32
X,Y=np.meshgrid(np.linspace(0.01*Lx,0.99*Lx,Nx),np.linspace(0.01*Ly,0.99*Ly,Ny))
u,v=VelocityField(pn,X,Y)
size=12
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
#plt.streamplot(X,Y,u,v,density=10,linewidth=1,arrowsize=1,arrowstyle='->')
plt.quiver(X,Y,u,v)
plt.plot(np.append([p.xa for p in pn],pn[0].xa),\
        np.append([p.ya for p in pn],pn[0].ya),\
        linestyle='-',linewidth=1,\
        marker='o',markersize=6,color='#CD2305');
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.title('Contour of velocity field');
plt.show()
