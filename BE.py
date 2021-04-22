# -*- coding: utf-8 -*-
"""
Created on Tue Apr 13 16:22:06 2021

@author: fbrav
"""

import numpy as np
from numpy import array
import matplotlib.pyplot as plt

def def_beam_cantilever(H,L,b,dens,E = 4e10):
    A = H*b
    I = b*H**3/12
    q = dens*A
    c1 = -q*L
    c2 = -q*L**2/2-c1*L
    c3 = 0
    c4 = 0
    delta = 1/(E*I)*(q*L**4/24+c1*L**3/6+c2*L**2/2+c3*L+c4)  
    return delta

def curv_beam_cantilever(H,L,b,dens,E = 4e10):
    A = H*b
    I = b*H**3/12
    q = dens*A
    c1 = -q*L
    c2 = -q*L**2/2-c1*L
    c3 = 0
    c4 = 0
    curv = 1/(E*I)*(q*0**2/2+c1*0+c2)  
    return curv


b = 1 #m
rho = 2400 #kg/m3
g = -9.81 #m/s2
fc= 250
E = 4700*np.sqrt(fc)*10**6
dens = rho*g

L0 = 19 #m
H0 = 1  #m
d0 = L0/H0
print(f"L0/H0 = {L0/H0}")
L1 = 19 #m
H1 = 1.2  #m
d1 = L1/H1
print(f"L1/H1 = {L1/H1}")
L2 = 19 #m
H2 = 1.8  #m
d2 = L2/H2
print(f"L2/H2 = {L2/H2}")
L3 = 19 #m
H3 = 2.5  #m
d3 = L3/H3
print(f"L3/H3 = {L3/H3}")
L4 = 19 #m
H4 = 4  #m
d4 = L4/H4
print(f"L4/H4 = {L4/H4}")
#L5 = 19 #m
#H5 = 6  #m
#d5 = L5/H5
#print(f"L5/H5 = {L5/H5}")


delta0 = def_beam_cantilever(H0,L0,b,dens,E)
delta1 = def_beam_cantilever(H1,L1,b,dens,E)
delta2 = def_beam_cantilever(H2,L2,b,dens,E)
delta3 = def_beam_cantilever(H3,L3,b,dens,E)
delta4 = def_beam_cantilever(H4,L4,b,dens,E)
#delta5 = def_beam_cantilever(H5,L5,b,dens,E)

sigma0 = -E*curv_beam_cantilever(H0,L0,b,dens,E)
sigma1 = -E*curv_beam_cantilever(H1,L1,b,dens,E)
sigma2 = -E*curv_beam_cantilever(H2,L2,b,dens,E)
sigma3 = -E*curv_beam_cantilever(H3,L3,b,dens,E)
sigma4 = -E*curv_beam_cantilever(H4,L4,b,dens,E)
#sigma5 = -E*curv_beam_cantilever(H5,L5,b,dens,E)

print(f"delta0 = {delta0} m")
print(f"delta1 = {delta1} m")
print(f"delta2 = {delta2} m")
print(f"delta3 = {delta3} m")
print(f"delta4 = {delta4} m")
#print(f"delta5 = {delta5} m")


print(f"sigma0 = {-E*curv_beam_cantilever(H0,L0,b,dens,E)}")
print(f"sigma1 = {-E*curv_beam_cantilever(H1,L1,b,dens,E)}")
print(f"sigma2 = {-E*curv_beam_cantilever(H2,L2,b,dens,E)}")
print(f"sigma3 = {-E*curv_beam_cantilever(H3,L3,b,dens,E)}")
print(f"sigma4 = {-E*curv_beam_cantilever(H4,L4,b,dens,E)}")
#print(f"sigma5 = {-E*curv_beam_cantilever(H5,L5,b,dens,E)}")



plt.figure()
plt.plot([d0,d1,d2,d3,d4], [delta0,delta1,delta2,delta3,delta4], 'ro')
plt.title('Deformación vs L/H')
plt.xlabel('(L/H)')
plt.ylabel('δ(L/H)')
plt.axis([0, 20, 0, -0.07])
plt.show()


#plt.figure(2)
plt.plot([d0,d1,d2,d3,d4], [sigma0,sigma1,sigma2,sigma3,sigma4], 'ro')
plt.title('Tesión vs L/H')
plt.xlabel('(L/H)')
plt.ylabel('σ(L/H)')
plt.axis([0, 20, 0, 54000000])
plt.show()
