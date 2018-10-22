import matplotlib.pyplot as plt
import scipy.integrate as integrate
import numpy as np


def analytic_potential(r0 = 1e-5):
    def coth(x):
        return 1/np.tanh(x)

    def sech(x):
        return 1/np.cosh(x)



    def phi0(xi):
        return 2 * np.log(np.abs(np.tanh(0.5 * (xi-xi_0))))


    def G(xi):
        return np.log(np.abs(np.cosh(0.5 * (xi-xi_0))/np.cosh(0.5 * (xi_0))))

    def I_G(xi):
        I = integrate.quad(lambda x: G(x), 0, xi)
        return I[0]

    def phiprime(xi):
        return xi ** 2 / 2 - 2 * g * xi + 2 * (2 * g - xi ) * coth(0.5*(xi-xi_0))

    def phiprime0(xi):
        return coth( 0.5 * (xi-xi_0) ) * ( sech( 0.5 * (xi-xi_0) ) ) ** 2

    def phi1(xi):
        aRes = []
        c = -phiprime(length)
        for i in range(0, int(len(xi))):
            xi1 = xi[i]
            res  = xi1 ** 3 / 6 - g * xi1 ** 2 + c * xi1 + 8 * g * (G(xi1) - G(0)) - 4 * G(xi1) + 4 * I_G(xi1)
            aRes.append(res)
        return np.array(aRes)

    def phi(xi):
        print("r value = " + str(r))
        return phi0(xi) + r/k * phi1(xi)

    def Efield(xi):
        return phiprime0(xi) + r * phiprime(xi)
    def Efield0(xi):
        return phiprime0(xi)

    def plot_potential(xi, p0, p1):
        plt.figure()
        plt.title('Potential to first order in r', fontsize=16, fontweight='bold')
        plt.xlabel(r'$\xi = \kappa x$', fontsize=16)
        plt.ylabel(r'Dimentionless potential ', fontsize=16)
        plt.plot(xi, p0, 'b', ms = 1, label='Zero order potential')
        plt.plot(xi, p1, 'g--', ms = 1, label='First order potential')

        plt.legend(loc='upper right', shadow=True, fontsize='x-large').get_frame().set_facecolor('#00FFCC')
        plt.subplots_adjust(hspace=0.4)

        plt.grid(True, color= '#F2F2F2')

        ## THIS COMMAND PRINTS THE FIGURE (COMPLETE) TO EPS
        plt.savefig('potentials.eps', format='eps', dpi=1000, fontsize=16, fontweight='bold')
        ## This command shows the plot.
        plt.show()



    #########################  PLOTTING SCRIPT  ###############################################
    #plt.style.use('seaborn')

    N = 1E4
    C_b = 0.1


    z = 2
    e = 1.60217662E-19
    k = 1.38064852E-23
    T = 300
    Na = 6.02E23
    Fa = Na * e
    R = Na * k
    coef = z*e/(k*T)
    epsilon = 80.9 * 8.85418782E-12
    k = np.sqrt(0.1 * (z*Fa) ** 2/ (R * T * epsilon))
    d = 1/k
    r = - k * r0
    V0 = -coef*0.15


    xi_0 = np.log(np.abs(np.tanh(V0/4)))
    g = coth(xi_0/2)
    length = 20 * k * d
    step = np.float_(length/N)
    xi = np.arange(0, length, step)

    p1 = phi(xi)
    p0 = phi0(xi)

    E1 = Efield(xi)
    E0 = Efield0(xi)

    #plot_potential(xi, p0, p1)

    file = open('../results/potential0.txt', 'w')
    file2 = open('../results/E-ana-r' + str(r0) + '.txt', 'w')
    file3 = open('../results/E1-ana-r' + str(r0) + '.txt', 'w')

    for i in range(0, len(xi)):
        file.write(str(xi[i]) + "\t" +  str(p0[i]) + "\n")
        file2.write(str(xi[i]) + "\t" + str(E0[i]) + "\n")
        file3.write(str(xi[i]) + "\t" + str(E1[i]) + "\n")

    file.close()

    file = open('../results/potential1-A'+str(r0)+'.txt', 'w')
    for i in range(0, len(xi)):
        file.write(str(xi[i]) + "\t" +  str(p1[i]) + "\n")

    file.close()
    file2.close()
    file3.close()


#reactRate = [5E-3, 3E-3, 5E-3]
reactRate = [5E-8, 3E-8, 1E-8]
for i in range(0,len(reactRate)):
    analytic_potential(reactRate[i])
