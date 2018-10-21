from numericMethods.ode.run_kut4 import integrate
from numericMethods.rootFinding.ridders import ridder
import numpy as np
from math import isnan



#Define the initial conditions as a function of u - in the shooting method we don't know the proper initial condition for E(0)
#such that the condition phi(d)=V0, where d is the laminar flux region.


#reactRate = [5E-3, 3E-3, 5E-3]
reactRate = [5E-8, 3E-8, 1E-8]
for p in range(0,len(reactRate)):

    def init(u):
        return np.array([0.1, 0.1, u, 0])

    # Define the function F(x,y) such that y' = F(x,y)
    def F(x,y):
        #C plus
        Cp =  y[0]
        #C minus
        Cm = y[1]
        #E electric field
        E = y[2]
        #Electric potential
        V = y[3]
        #Function
        #r0 = r/k
        F = np.zeros(4)
        F[0] = Cp*E - r0
        F[1] = -Cm*E
        F[2] = (Cp-Cm)
        F[3] = -E
        return F

    def invert(A):
        B = []
        N = len(A)
        for i in range(0,N):
            B.append(A[N-1-i])
        return np.array(B)


    def find_not_nan(L):
        for i in range(0,len(L)-1):
            if(isnan(L[len(L)-1-i]) == False):
                return len(L)-1-i

    def r(u):
        X,Y = integrate(F, xStart, init(u), xStop, step)
        V_u = Y[-1,3]
        if(isnan(V_u)):
            index = int(find_not_nan(Y[:,3]))
            V_u = Y[index,3]

        return (V_u-V_0)/V_0

    #Here we start the algorithm to find de numeric solution to the system defined by y'=F(x,y). Note that here y
    #is a vector with four components

    #Define integration range.

    xStop = 20.0
    xStart = 0.0

    #Define the number of steps of integration.
    N = 1E3
    step = (xStop-xStart)/N
    r0 = reactRate[p]
    u1 = 10
    u2 = 0

    z = 2
    e = 1.60217662E-19
    k = 1.38064852E-23
    T = 300
    Na = 6.02E23
    Fa = Na * e
    R = Na * k
    coef = z*e/(k*T)
    V_0 = -coef*0.15

    epsilon = 80.9 * 8.85418782E-12
    k = np.sqrt(0.1 * (z*Fa) ** 2/ (R * T * epsilon))

    u = ridder(r, u1, u2)

    #Plot results, if initial condition for E(0) is found.


    import matplotlib.pyplot as plt

    if (u != None):
        X,Y = integrate(F, xStart, init(u), xStop, step)
        Cp = invert(Y[:,0])
        Cm = invert(Y[:,1])
        E = invert(Y[:,2])
        phi = invert(Y[:,3])
        xi = X
        print("Border condition at interface is: " + str(-77.357*0.15))
        print("Obtained border condition with shooting method: " + str(phi[0]))
        print("Obtained border condition in the bulk: " + str(E[-1]))
        print("Shooting error: " + str((phi[0]-V_0)/V_0))
        

    #Plot electric potential
        plt.figure(1)
        plt.title('Numeric Potential', fontsize=16, fontweight='bold')
        plt.xlabel(r'$\xi = \kappa x$', fontsize=16)
        plt.ylabel(r'Dimentionless potential ', fontsize=16)
        plt.plot(X[:], phi, 'g', label = "Potential")
        plt.axis([xStart, xStop, -20, -0.0])
        plt.legend()
        plt.savefig('potential.eps', format='eps', dpi=1000, fontsize=16, fontweight='bold')
        #plt.show()
        
        #Plot electric field
        plt.figure(2)
        plt.plot(X[:], E, 'g', label = "Electric Field")
        plt.axis([xStart, xStop, 0.0, 100])
        plt.legend()
        plt.savefig('Efield.eps', format='eps', dpi=1000, fontsize=16, fontweight='bold')
        #plt.show()

        #Plot Concetration
        plt.figure(3)
        plt.title('Numeric Concentration', fontsize=16, fontweight='bold')
        plt.xlabel(r'$\xi = \kappa x$', fontsize=16)
        plt.ylabel(r'Molar Concentration ', fontsize=16)
        plt.plot(X[:], Cp, 'b', label = "Concentration +")
        plt.axis([xStart, xStop, 0, 1])
        plt.plot(X[:], Cm, 'r', label = "Concentration -")
        plt.legend()
        
        #plt.savefig('Concentration-num-'+str(r0)+'.eps', format='eps', dpi=1000, fontsize=16, fontweight='bold')

        #plt.show()

        file = open('./results/potential-num-r'+str(r0)+'.txt', 'w')
        file2 = open('./results/cp-num-r'+str(r0)+'.txt', 'w')
        file3 = open('./results/cm-num-r'+str(r0)+'.txt', 'w')
        file4 = open('./results/E-num-r' + str(r0) + '.txt', 'w')


        for i in range(0, len(xi)):
            file.write(str(xi[i]) + "\t" +  str(phi[i]) + "\n")
            file2.write(str(xi[i]) + "\t" +  str(Cp[i]) + "\n")
            file3.write(str(xi[i]) + "\t" +  str(Cm[i]) + "\n")
            file4.write(str(xi[i]) + "\t" + str(E[i]) + "\n")
        file.close()
        file2.close()
        file3.close()
        file4.close()


