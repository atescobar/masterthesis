import numpy as np
from numpy.linalg import inv
import scipy as sp
import scipy.sparse

def get_static_sol(filename):
    file = open(filename, 'r')
    X = []
    Y = []
    for line in file:
        aux1, aux2 = line.split()
        X.append(float(aux1))
        Y.append(float(aux2))
    return X, Y

def get_next_Psi(Psi, A, b):
    invA = inv(A)
    return np.dot(invA, b)

def invert(A):
    aux = []
    for i in range(0,len(A)):
        aux.append(A[len(A)-1-i])
    return aux


#Adquirir condiciones de borde
print("retrieving files")
x_cp, numCp = get_static_sol('./borderConditions/cp-num-r5e-08.txt')
x_cm, numCm = get_static_sol('./borderConditions/cm-num-r5e-08.txt')
x_psi, numPsi = get_static_sol('./borderConditions/potential-num-r5e-08.txt')


Psi0 = numPsi[0]
dx = 1
Cb = 0.1

#Inicializar los valores
N = 500
M = len(x_cp)

Cp = np.zeros([N,M])
Cm = np.zeros([N,M])
Psi  = np.zeros([N,M])
print("setting up border condition")
for j in range(1,M):
    Cp[0, j] = numCp[j]
    Cm[0, j] = numCm[j]
    Psi[0, j] = numPsi[j]
    #Cp[0, j] = Cb
    #Cm[0, j] = Cb
    #Psi[0, j] = 0

for i in range(0,N):
    Cp[i,0] = Cb
    Cm[i,0] = Cb
    Psi[i,0] = Psi0
    Psi[i,M-1] = 0


#hasta aca asumo que conocemos las condiciones de borde.
#definicion de constantes
z = 2
e = 1.60217662E-19
k = 1.38064852E-23
T = 300
Na = 6.02E23
Fa = Na * e
R = Na * k
coef = z * e / (k * T)
V_0 = -coef * 0.15

epsilon = 80.9 * 8.85418782E-12
k = np.sqrt(0.1 * (z * Fa) ** 2 / (R * T * epsilon))
dt = 1
k = dx * k
Dp = 0.00005
Dm = 0.00005
rp = dt * Dp / dx ** 2
rm = dt * Dm / dx ** 2

#r = k * 5E-8
r = 5E-8

#Aca definimos el sistema que resuleve para cada iteracion del potencial.
# La matriz de coeficientes es M x M donde M es el numero de pasos en el eje x
print("creating coefficient matrix")
a = np.ones(M)
b1 = -2 * np.ones(M)
c = np.ones(M)
positions = [-1, 0, 1]
A = sp.sparse.spdiags(np.array([a, b1, c]), positions, M, M).todense()
Ainv = np.linalg.inv(A)


print("done creating coefficient matrix")
print("Setting up resulting vector")
b = np.zeros([N,M])

print("done...")
print("System is ready to be solved... starting iteration now")
for n in range(0,N-1):
    for k in range(1, M-1):
        #compute next value of concentrations using previous ones
        Cp[n + 1, k] = rp * Cp[n, k + 1] * (1 + Psi[n, k] - Psi[n, k + 1]) + Cp[n, k] * (1 + rp * (-2 + Psi[n, k] - Psi[n, k - 1])) + rp * Cp[n, k - 1]
        Cm[n + 1, k] = rm * Cm[n, k + 1] * (1 + Psi[n, k + 1] - Psi[n, k]) + Cm[n, k] * (1 + rm * (-2 + Psi[n, k - 1] - Psi[n, k])) + rm * Cm[n, k - 1]
        #Concentration border condition
        if (k == M - 1):
            Cp[n + 1, k] = Cp[n + 1, k] - dt * r

        #Psi border condition
        if(k == 1 ):
            b[n+1,k] = -k ** 2 * (Cp[n+1,k] - Cm[n+1, k]) - Psi[n+1, 0]
        else:
            b[n+1,k] = -k ** 2 * (Cp[n+1,k] - Cm[n+1, k])

    print("\rCompleted: " + "time axis: " +str("{0:.2f}".format((n / N) * 100)) + " %", end="")
    Psi[n+1,:] = np.dot(Ainv, b[n+1, :])

print("\n done computing")
#print("writing results files...")
### Write files
#file1 = open('dynamicResults/result1.txt', 'w')
#file2 = open('dynamicResults/result2.txt', 'w')
#file3 = open('dynamicResults/result3.txt', 'w')
#file4 = open('dynamicResults/result4.txt', 'w')

#for i in range(1, M):
#    file1.write(str(x_psi[i]) + '\t' + str(Psi[20,i]))
#    file2.write(str(x_psi[i]) + '\t' + str(Psi[40, i]))
#    file3.write(str(x_psi[i]) + '\t' + str(Psi[60, i]))
#    file4.write(str(x_psi[i]) + '\t' + str(Psi[80, i]))

#file1.close()
#file2.close()
#file3.close()
#file4.close()

#print("done \n Results stored in ./dynamicalResults" )


####################################################################################################
import matplotlib.pyplot as plt
plt.xlabel("dimensionless length paramenter")
plt.ylabel("dimensionless potentialr")
plt.plot(x_psi, Psi[0, :], label='Phi at 0')
#plt.plot(x_psi, Psi[1,:], label='Phi at 1')
#plt.plot(x_psi, Psi[3,:], label='Phi at 2')
plt.show()
