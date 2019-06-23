import math
from numpy import sign
from numpy import abs
import numpy as np
def ridder(f,a,b,tol=1e-2):
    a = np.array(a)
    b = np.array(b)
    fa = f(a)
    fb = f(b)
    if ((sign(fb) == sign(fa)).all()):
        print('Root not bracketed :(')
        return None

    for i in range(400):
      	# Compute the improved root x from Ridder’s formula
        c = 0.5*(a + b)
        fc = f(c)
        s = np.sqrt(fc**2 - fa*fb)
        if s.any() == 0.0: 
            return None
        print('fa: '+str(f(a)))
        print('fb: '+str(f(b)))
        print('fc: '+str(f(c)))
        dx = (c - a)*fc/s
        for k in range(len(fa)):
            if (fa[k] - fb[k]) < 0.0: 
                dx = -dx

        x = c + dx
        fx = f(x)
        
      	# Test for convergence: last increment must be less than the predefined tolerance in units of the 
      	# last known root value x.
        if i > 0:
            if (abs(fx).any() < tol):#(x - xOld < tol*max(abs(x),1.0)):
                print('Root found!: ' + 'x = ' +  str(xOld) + " => " + 'r(x) = ' + str(fx))
                print('Tolerance: ' + str(tol) )
                return x

        xOld = x

      	# Re-bracket the root as tightly as possible
        for k in range(len(fc)):
            if sign(fc[k]) == sign(fx[k]):
                if sign(fa[k])!= sign(fx[k]): 
                    b[k] = x[k]
                    fb[k] = fx[k]

                else: 
                    a[k] = x[k] 
                    fa[k] = fx[k]
            else:
                a[k] = c[k]
                b[k] = x[k]
                fa[k] = fc[k]
                fb[k] = fx[k]

    return None
    print('Too many iterations')
