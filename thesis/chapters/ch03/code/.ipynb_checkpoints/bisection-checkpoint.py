import math
from numpy import sign
from numpy import abs

def findRoot(f, m, a, b, tol=1e-5):
    if (sign(f(a, m)) == sign(f(b,m))):
        print("Root not bracketed")
        return
    if(a > b):
#        print("a=" + str(a) +" was bigger than b=" + str(b) +". Will switch them")
        aux = b
        b = a
        a = aux

    
    intervalLength = b - a
    i=0
    while(abs(intervalLength) > tol):
        i += 1
        c = (b + a)/2
        if ( sign(f(c,m)) == sign(f(a,m)) ):
            a = c
        if ( sign(f(c,m)) == sign(f(b,m)) ):
            b = c
        intervalLength = b - a

        if (abs(intervalLength) < tol):
            return c


        
        
    