from READING import sigma, desp
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



delta0 = def_beam_cantilever(H0,L0,b,dens,E)
delta1 = def_beam_cantilever(H1,L1,b,dens,E)
delta2 = def_beam_cantilever(H2,L2,b,dens,E)
delta3 = def_beam_cantilever(H3,L3,b,dens,E)
delta4 = def_beam_cantilever(H4,L4,b,dens,E)


sigma0 = -E*curv_beam_cantilever(H0,L0,b,dens,E)
sigma1 = -E*curv_beam_cantilever(H1,L1,b,dens,E)
sigma2 = -E*curv_beam_cantilever(H2,L2,b,dens,E)
sigma3 = -E*curv_beam_cantilever(H3,L3,b,dens,E)
sigma4 = -E*curv_beam_cantilever(H4,L4,b,dens,E)



s1 = sigma0/sigma("sigma_x10-025.msh")
s2 = sigma0/sigma("sigma_x10-04.msh")
s3 = sigma0/sigma("sigma_x10-10.msh")

s4 = sigma1/sigma("sigma_x12-025.msh")
s5 = sigma1/sigma("sigma_x12-04.msh")
s6 = sigma1/sigma("sigma_x12-10.msh")

s7 = sigma2/sigma("sigma_x18-025.msh")
s8 = sigma2/sigma("sigma_x18-04.msh")
s9 = sigma2/sigma("sigma_x18-10.msh")

s10 = sigma3/sigma("sigma_x25-025.msh")
s11 = sigma3/sigma("sigma_x25-04.msh")
s12 = sigma3/sigma("sigma_x25-10.msh")

s13 = sigma4/sigma("sigma_x40-025.msh")
s14 = sigma4/sigma("sigma_x40-04.msh")
s15 = sigma4/sigma("sigma_x40-10.msh")

d1 = delta0/desp("desplazamientos10-025.msh")
d2 = delta0/desp("desplazamientos10-04.msh")
d3 = delta0/desp("desplazamientos10-10.msh")

d4 = delta1/desp("desplazamientos12-025.msh")
d5 = delta1/desp("desplazamientos12-04.msh")
d6 = delta1/desp("desplazamientos12-10.msh")

d7 = delta2/desp("desplazamientos18-025.msh")
d8 = delta2/desp("desplazamientos18-04.msh")
d9 = delta2/desp("desplazamientos18-10.msh")

d10 = delta3/desp("desplazamientos25-025.msh")
d11 = delta3/desp("desplazamientos25-04.msh")
d12 = delta3/desp("desplazamientos25-10.msh")

d13 = delta4/desp("desplazamientos40-025.msh")
d14 = delta4/desp("desplazamientos40-04.msh")
d15 = delta4/desp("desplazamientos40-10.msh")
    



plt.plot([d1,d2,d3], [s1,s2,s3],"ro", )
plt.plot([d4,d5,d6], [s4,s5,s6],"bo")
plt.plot([d7,d8,d9], [s7,s8,s9],"ko")
plt.plot([d10,d11,d12], [s10,s11,s12],"mo")
plt.plot([d13,d14,d15], [s13,s14,s15],"yo")

plt.title("FEM - BE")
plt.xlabel('σ_FEM(L/H)/σ_BE(L/H)')
plt.ylabel('δ_BE(L/H)/δ_FEM(L/H)')
plt.axis([0, 50, 0, 4.0])
plt.legend(["L0/H0 = 19.0","L1/H1 = 15.83","L2/H2 = 10.55","L3/H3 = 7.6","L4/H4 = 4.75"])
plt.savefig("graficoBE-FEM.png")
plt.show()

#plt.plot([d0,d1,d2,d3,d4,d5], [sigma0,sigma1,sigma2,sigma3,sigma4, sigma5], 'ro')
#plt.xlabel('(L/H)')
#plt.ylabel('σ(L/H)')
#plt.axis([0, 20, 0, 51000000])
#plt.show()

plt.plot([d1,d4,d7,d10,d13], [s1,s4,s7,s10,s13],"1")
plt.plot([d2,d5,d8,d11,d14], [s2,s5,s8,s11,s14],"x")
plt.plot([d3,d6,d9,d12,d15], [s3,s6,s9,s12,s15],"*")
plt.title("FEM - BE")
plt.xlabel('σ_FEM(L/H)/σ_BE(L/H)')
plt.ylabel('δ_BE(L/H)/δ_FEM(L/H)')
plt.axis([0, 50, 0, 4.0])
plt.legend(["b=0.25","b=0.4","b=1.0"])
plt.savefig("grafico grilla BE-FEM.png")
plt.show()
