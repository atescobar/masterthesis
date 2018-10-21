from run_kut4 import integrate
from falsePosition import falsePosition
import numpy as np

# We want to solve a system of the type y' = F(x,y)


#Define function F

def F(x,y):
	F = np.zeros(4)
	a = 0.1
	cp = y[0]
	cm = y[1]
	E = y[2]
	phi = y[3]
	
	if(E == None):
		print("function argument error: value of E is None")
	if(phi == None):
		print("function argument error: value of phi is None")
	if(cp == None):
		print("function argument error: value of cpis None")
	if(cm == None):
		print("function argument error: value of cm is None")
		
	F[0] = 1.2851E-22*cp*E-0.1
	F[1] = -1.2851E-22*cm*E
	F[2] = 4.5176E-10*(cp-cm)
	F[3] = -E

	return F

def r(u):
	print("calling residual function")
	N = 1E4
	step = (xStop-xStart)/N
	X, Y = integrate(F, xStart, init(u), xStop, step)
	if (Y[-1,3]):
		print( "value of r :=" + str(Y[-1,3] + 0.15) )
		res = Y[-1,3] + 0.15
		return res
	else:
		print("runge-kutta method failed to integrate")
	
def init(u):
	Na = 6.02E23
	return np.array([0.1*Na, 0.1*Na, u, 0.0])

xStop = 4.2603596177667828e-06
xStart = 0.0
N = 1E4

step = (xStop-xStart)/N
u1 = 10
u2 = -10
print(r(u1))
print(r(u2))
if (r(u1)*r(u2)<0):
	print("wena choro")
u = falsePosition(r, u1, u2)
X,Y = integrate(F, xStart, init(u), xStop, step)


import matplotlib.pyplot as plt

print("Y0\tY1\tY2\tY3")

for i in range(0,20):
	print(Y[0,:])
	
	
plt.plot(X[:], Y[:,1], 'r')
plt.axis([0.0, 4.2603596177667828e-06, -0.2, -0.0])
plt.show()
