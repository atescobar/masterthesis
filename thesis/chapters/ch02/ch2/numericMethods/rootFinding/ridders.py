import math
from numpy import sign
from numpy import abs

def ridder(f,a,b,tol=1e-2):
    fa = f(a)
    if fa == 0.0: 
    	return a


    fb = f(b)
    if fb == 0.0: 
    	return b


    if (sign(fb) == sign(fa)):
        print('Root not bracketed :(')
        return None

    for i in range(400):
      	# Compute the improved root x from Ridder’s formula
        c = 0.5*(a + b)
        fc = f(c)
        s = math.sqrt(fc**2 - fa*fb)
        if s == 0.0: 
            return None

        dx = (c - a)*fc/s

        if (fa - fb) < 0.0: 
        	dx = -dx

        x = c + dx
        fx = f(x)

      	# Test for convergence: last increment must be less than the predefined tolerance in units of the 
      	# last known root value x.
        if i > 0:
            if (abs(fx) < tol):#(x - xOld < tol*max(abs(x),1.0)):
                print('Root found!: ' + 'x = ' +  str(xOld) + " => " + 'r(x) = ' + str(fx))
                print('Tolerance: ' + str(tol) )
                return x

        xOld = x

      	# Re-bracket the root as tightly as possible
        if sign(fc) == sign(fx):
            if sign(fa)!= sign(fx): 
            	b = x
            	fb = fx

            else: 
            	a = x 
            	fa = fx
        else:
            a = c
            b = x
            fa = fc
            fb = fx

    return None
    print('Too many iterations')
