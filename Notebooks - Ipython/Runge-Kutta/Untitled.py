from run_kut4 import integrate
from falsePosition import falsePosition
import numpy as np

# We want to solve a system of the type y' = F(x,y)


#Define function F
#y'=F(x,y)
#y=(y',y)=(y1,y2)

def F(x,y):
	F = np.zeros(2)
	F[0] = -4*y[1]
	F[1] = y[0]
	return F

#residual: y'(0) such that y(pi/2)=0.0?
def r(u):
	N = 100
	step = (xStop-xStart)/N
	X, Y = integrate(F, xStart, init(u), xStop, step)
	res = Y[-1,0]
	return res

def init(u):
	return np.array([u, 1.0])

xStop = np.pi/2
xStart = 0.0
N = 100
step = (xStop-xStart)/N

u1 = 5
u2 = -5
print(r(0))

u = falsePosition(r, u1, u2)
X,Y = integrate(F, xStart, init(u), xStop, step)


import matplotlib.pyplot as plt
	
	
plt.plot(X[:], Y[:,1], 'r')
plt.axis([0.0, np.pi/2, -1.0, 1.0])
plt.show()
