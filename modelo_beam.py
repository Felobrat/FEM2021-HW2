import numpy as np
from sympy import *
import matplotlib.pylab as plt
from cst import cst_planestress, cst_planestress_post

fid = open("beam.msh", "r")

LINE_ELEMENT = 1
TRI_ELEMENT = 2
FIJO = 1
VIGA = 2

while True:
	line = fid.readline()

	if line.find("Nodes") >= 0:
		break
Nnodes = int(fid.readline())

xy = np.zeros([Nnodes, 2])

for i in range(Nnodes):
	line = fid.readline()
	sl = line.split()
	xy[i, 0] = float(sl[1])
	xy[i, 1] = float(sl[2])

print(xy) 

while True:
	line = fid.readline()

	if line.find("$Elements") >= 0:
		break
Nelements = int(fid.readline())

print(Nelements)

conec = np.zeros([Nelements,3], dtype = np.int32)

fixed_nodes = []
Ntriangles = 0
Triangless = []

for i in range(Nelements):
	line = fid.readline()
	sl = line.split()
	element_number = np.int32(sl[0]) -1 
	element_type = np.int32(sl[1])
	physical_grp = np.int32(sl[3])
	entity_number = np.int32(sl[4])

	if element_type == LINE_ELEMENT and \
	 physical_grp == LINE_ELEMENT: 
		n1 = np.int32(sl[5]) -1
		n2 = np.int32(sl[6]) -1
		fixed_nodes += [n1,n2]

	if element_type == TRI_ELEMENT and \
	 physical_grp == TRI_ELEMENT:
		 n0 = np.int32(sl[5]) -1
		 n1 = np.int32(sl[6]) -1
		 n2 = np.int32(sl[7]) -1 
		 
		 conec[element_number, :] = [n0, n1, n2]

		 Triangless.append(element_number)
		 Ntriangles +=1

print(conec)

NDOFs = 2*Nnodes 

K = np.zeros((NDOFs, NDOFs))
f = np.zeros((NDOFs, 1))

properties = {}
rho = 2500. 
g = 9.81	
fc= 250

properties["E"] = 4700*np.sqrt(fc)*10**6
properties["v"] = 0.25
properties["bx"] = 0
properties["by"] = -rho*g
properties["t"] = 1. 


for e in Triangless:  #ingresa cada elemento e con su respecitov ke, fe y lo ensambla con el DSM
	ni = conec[e,0]
	nj = conec[e,1]
	nk = conec[e,2]

	print(f"e = {e} ni = {ni} nj = {nj} nk = {nk}")

	xy_e = xy[[ni, nj, nk],:]

	#print(f"xy_e = {xy_e}")

	ke, fe = cst_planestress(xy_e, properties)

	#print(f"ke = {ke}")

	#Node k --> [3*k, 3*k+1, 3*k+2]

	d = [2*ni, 2*ni+1, 2*nj, 2*nj+1, 2*nk, 2*nk+1]  #esto es lo que cambia con otro tipo de elemento, desplazamiento hori, ver, giro por nodo 

	#DSM
	for i in range(len(d)): #limite de loops también cambian con otro elemento 
		p = d[i]
		for j in range(len(d)):
			q = d[j]
			K[p,q] += ke[i,j]
		f[p] += fe[i]

#hasta aqui tenemos emsablada la matriz de rigidez 

fixed_nodes = np.unique(fixed_nodes)

constrained_DOFs = []

for n in fixed_nodes:
	constrained_DOFs += [2*n, 2*n+1]

free_DOFs = np.arange(NDOFs)
free_DOFs = np.setdiff1d(free_DOFs,constrained_DOFs )

print(f"fixed_nodes = {fixed_nodes}")
print(f"constrained_DOFs = {constrained_DOFs}")
print(f"free_DOFs = {free_DOFs}")

Kff = K[np.ix_(free_DOFs, free_DOFs)]
Kfc = K[np.ix_(free_DOFs, constrained_DOFs)]
Kcf = K[np.ix_(constrained_DOFs, free_DOFs)]
Kcc = K[np.ix_(constrained_DOFs, constrained_DOFs)]

ff = f[free_DOFs]
fc = f[constrained_DOFs]

# Solve
from scipy.linalg import solve
u = np.zeros((NDOFs,1))

u[free_DOFs] = solve(Kff, ff)

#Get reaction forces 
R = Kcf @ u[free_DOFs] + Kcc @ u[constrained_DOFs]	- fc 

print (f"u = {u}")
print (f"R = {R}")


factor = 1e2

uv= u.reshape([-1, 2])  #*factor

plt.plot(xy[:,0] + factor*uv[:,0], xy[:,1]+ factor*uv[:,1], '.')

for e in Triangless:  #ingresa cada elemento e con su respecitov ke, fe y lo ensambla con el DSM
	ni = conec[e,0]
	nj = conec[e,1]
	nk = conec[e,2]
	xy_e = xy[[ni, nj, nk, ni],:] + factor*uv[[ni, nj, nk, ni],:]  #esto para muchos elementos se va a trabar en matplotlib, visualizar en gmsh

	plt.plot(xy_e[:,0], xy_e[:,1], 'k')

plt.axis('equal')
plt.show()

#
plt.matshow(K)
plt.show()

#AQUI EMPIEZA EL POST PROCESO 

from gmsh_post import write_node_data, write_node_data_2, write_element_data

nodes = np.arange(1, Nnodes+1)

write_node_data("ux.msh", nodes, uv[:,0], "Despl. X" )
write_node_data("uy.msh", nodes, uv[:,1], "Despl. Y" )

write_node_data_2("desplazamientos.msh", nodes, uv[:,0] ,uv[:,1], "Despl" )

#calculo de tensiones

sigmaxx = np.zeros(Ntriangles+1)
sigmayy = np.zeros(Ntriangles+1)
sigmaxy = np.zeros(Ntriangles+1)

i = 0
for e in Triangless:  #ingresa cada elemento e con su respecitov ke, fe y lo ensambla con el DSM
	ni = conec[e,0]
	nj = conec[e,1]
	nk = conec[e,2]
	xy_e = xy[[ni, nj, nk, ni],:] 

	uv_e = uv[[ni, nj, nk],:]  #esto para muchos elementos se va a trabar en matplotlib, visualizar en gmsh

	u_e = uv_e.reshape((-1)) 

	Ee,  sigmae = cst_planestress_post(xy_e, u_e, properties)

	sigmaxx[i] = sigmae[0] 
	sigmayy[i] = sigmae[1]
	sigmaxy[i] = sigmae[2]

	i += 1

	#plt.plot(xy_e[:,0], xy_e[:,1], 'k')


elements = np.array(Triangless)+1
write_element_data("sigma_x.msh", elements, sigmaxx, "Sigma_x" )
write_element_data("sigma_y.msh", elements, sigmayy, "Sigma_y" )
write_element_data("sigma_xy.msh", elements, sigmaxy, "Sigma_xy" )





