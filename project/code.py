import numpy as np
from scipy import integrate
from math import *
import matplotlib.pyplot as plt

Ni = 16;       # Number of elements per side
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

class Panel:
    def __init__(self,xa,ya,xb,yb,ielm):
        self.xa,self.ya = xa,ya                     # 1st end-point
        self.xb,self.yb = xb,yb                     # 2nd end-point
        self.xc,self.yc = (xa+xb)/2,(ya+yb)/2       # control point
        self.ielm=ielm              #0=horizontal element,1=vertical element
        # location of the panel
        if (self.beta<=pi): self.loc = 'extrados'
        else: self.loc = 'intrados'
