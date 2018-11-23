#imports
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags

# Define grid parameters
D = 1.07 # Diffusion Coefficient
N = 100
tau = np.linspace(0,1, N) #shape is N+1

dt = 1/(N)  # N Partitions
dx = 1/(M) # N Partitions 
a = dt / dx ** 2 * D

# Define the coefficient matrix
di = ( 1 + 2 * a ) * np.ones(M-2)
di[0] = ( 1 + a )
A = diags(np.array([- a * np.ones(M-3), di, - a * np.ones(M-3)]), [-1, 0, 1], shape=(M-2, M-2)).toarray()
A_inv = np.asarray(np.linalg.inv(A))



def C(t):
    # Set up initial conditions for \rho
    Cb = 100

    rho = np.zeros([N, M])

    rho[0, :] = - Cb
    rho[0, :] = - Cb
    rho[:, -1] = 0    
    rho[:, 0] = rho[:, 1]
    
    #Starting iteration
    for n in range(0, N-1):
        rho[n+1, 1:M-1] = np.matmul(A_inv, rho[n, 1:M-1])
        rho[n+1, 0] = rho[n+1, 1]
        
    n = int(t/dt)
    
    return Cb * np.ones(M) + rho[n, :]