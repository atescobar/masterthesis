### In this script we try to find the solution to the problem
##
## y'' + 4 y = 0
## y(0) = 1.0
## y(pi/2) = 0

from numericMethods.ode.run_kut4 import integrate
from numericMethods.rootFinding.ridders import ridder
import numpy as np

### We want to solve a system of the type y' = F(x,y)
## Define F
def F(x,y):
	F = np.zeros(2)
	F[0] = -4*y[1]
	F[1] = y[0]
	return F

## Define the residual: looking for y'(0) such that y(pi/2)=0.0.
def r(u):
	N = 100
	step = (xStop-xStart)/N
	X, Y = integrate(F, xStart, init(u), xStop, step)
	res = Y[-1,0]
	return res

def init(u):
	return np.array([u, 1.0])

## Define intervale of interest and step.
xStop = np.pi/2
xStart = 0.0
N = 100
step = (xStop-xStart)/N

## Define initial guess for the y'(0) condition. It hast to bracket the actual root of the residual function r(u).
u1 = 5
u2 = -5

## Use ridder method for finding the root (Newton Raphson is avoided since we would need to find r'(u), which can't be obtained analytically.
u = ridder(r, u1, u2)

## Once u is found, we integrate one last time with Runge Kuta method.
X,Y = integrate(F, xStart, init(u), xStop, step)


## Plot the results
import matplotlib.pyplot as plt
	
	
plt.plot(X[:], Y[:,1], 'r', label='Numeric Solution')
plt.plot(X[:], np.cos(2*X[:]), 'b--', label='Analytic Solution')
plt.axis([0.0, np.pi/2, -1.0, 1.0])
plt.legend()
plt.show()
