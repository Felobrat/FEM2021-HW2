import numpy as np
from numpy import array, sqrt


def cst_planestress(xy, properties):  #nos casamos con una de las formulaciones
	
	E = properties["E"]
	v = properties["v"]
	bx = properties["bx"]
	by = properties["by"]
	t = properties["t"]

	Eσ = E / (1-v**2) * np.array(
		[
		[1, v, 0],
		[v, 1, 0],
		[0, 0, (1-v)/2]
		]) 


	x0 = xy[0,0]
	x1 = xy[1,0]
	x2 = xy[2,0]
	y0 = xy[0,1]
	y1 = xy[1,1]
	y2 = xy[2,1]

	l0 = np.sqrt((x1- x2)**2 + (y1- y2)**2)
	l1 = np.sqrt((x0- x2)**2 + (y0- y2)**2)
	l2 = np.sqrt((x1- x0)**2 + (y1- y0)**2)

	# formula de heron de alejandria 
	s = (l0 +l1 +l2)/2
	Ae = sqrt( (s-l0)*(s-l1)*(s-l2)*s )

	dzeta0_dx = (y1 - y2)/(x0*y1 - x0*y2 - x1*y0 + x1*y2 + x2*y0 - x2*y1)
	dzeta0_dy = (-x1 + x2)/(x0*y1 - x0*y2 - x1*y0 + x1*y2 + x2*y0 - x2*y1)
	dzeta1_dx = (-y0 + y2)/(x0*y1 - x0*y2 - x1*y0 + x1*y2 + x2*y0 - x2*y1)
	dzeta1_dy = (x0 - x2)/(x0*y1 - x0*y2 - x1*y0 + x1*y2 + x2*y0 - x2*y1)
	dzeta2_dx = (y0 - y1)/(x0*y1 - x0*y2 - x1*y0 + x1*y2 + x2*y0 - x2*y1)
	dzeta2_dy = (-x0 + x1)/(x0*y1 - x0*y2 - x1*y0 + x1*y2 + x2*y0 - x2*y1)

	B = array([
	[dzeta0_dx, 0 , dzeta1_dx, 0 , dzeta2_dx, 0],
	[0, dzeta0_dy, 0 , dzeta1_dy, 0 , dzeta2_dy],
	[dzeta0_dy, dzeta0_dx, dzeta1_dy, dzeta1_dx, dzeta2_dy, dzeta2_dx]	
	])

	print(Ae)

	ke = B.T @ Eσ @ B * Ae * t 
	fe = (Ae * t / 3) * array([bx, by, bx, by, bx, by])


	return ke, fe 

def cst_planestress_post(xy, u_e, properties):  #nos casamos con una de las formulaciones
	
	E = properties["E"]
	v = properties["v"]
	bx = properties["bx"]
	by = properties["by"]
	t = properties["t"]

	Eσ = E / (1-v**2) * np.array(
		[
		[1, v, 0],
		[v, 1, 0],
		[0, 0, (1-v)/2]
		]) 


	x0 = xy[0,0]
	x1 = xy[1,0]
	x2 = xy[2,0]
	y0 = xy[0,1]
	y1 = xy[1,1]
	y2 = xy[2,1]

	l0 = np.sqrt((x1- x2)**2 + (y1- y2)**2)
	l1 = np.sqrt((x0- x2)**2 + (y0- y2)**2)
	l2 = np.sqrt((x1- x0)**2 + (y1- y0)**2)

	# formula de heron de alejandria 
	s = (l0 +l1 +l2)/2
	Ae = sqrt( (s-l0)*(s-l1)*(s-l2)*s )

	dzeta0_dx = (y1 - y2)/(x0*y1 - x0*y2 - x1*y0 + x1*y2 + x2*y0 - x2*y1)
	dzeta0_dy = (-x1 + x2)/(x0*y1 - x0*y2 - x1*y0 + x1*y2 + x2*y0 - x2*y1)
	dzeta1_dx = (-y0 + y2)/(x0*y1 - x0*y2 - x1*y0 + x1*y2 + x2*y0 - x2*y1)
	dzeta1_dy = (x0 - x2)/(x0*y1 - x0*y2 - x1*y0 + x1*y2 + x2*y0 - x2*y1)
	dzeta2_dx = (y0 - y1)/(x0*y1 - x0*y2 - x1*y0 + x1*y2 + x2*y0 - x2*y1)
	dzeta2_dy = (-x0 + x1)/(x0*y1 - x0*y2 - x1*y0 + x1*y2 + x2*y0 - x2*y1)

	B = array([
	[dzeta0_dx, 0 , dzeta1_dx, 0 , dzeta2_dx, 0],
	[0, dzeta0_dy, 0 , dzeta1_dy, 0 , dzeta2_dy],
	[dzeta0_dy, dzeta0_dx, dzeta1_dy, dzeta1_dx, dzeta2_dy, dzeta2_dx]	
	])


	epsilon = B @ u_e
	sigma = Eσ @ epsilon


	return epsilon, sigma 


xy = array([
	[0,0],
	[1,0],
	[0,1]
	])

properties = {}
properties["E"] = 1.0
properties["v"] = 0.25
properties["bx"] = 0
properties["by"] = 1.0
properties["t"] = 1.0

ke, fe = cst_planestress(xy, properties)

print(f"ke = {ke}")
print(f"fe = {fe}")




