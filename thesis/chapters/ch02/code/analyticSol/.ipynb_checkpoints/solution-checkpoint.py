import numpy as np

params0 = {
    "z": 2,
    "e": 1.60217662E-19,
    "kb": 1.38064852E-23,
    "T": 300,
    "Na": 6.02E23,
    "Fa": 96485.3329,#Na * e
    "R": 8.314472,
    "V0": -0.15,
    "D1": 1.05,
    "D2":  0.76,
    "Cb": 1,
    #d = 1.544E-6
    "epsilon": 80.9 * 8.85418782E-12,
    "length": 1.0
}

class model:

    def __init__(self, xi, params = params0):
        self.params = params
        self.z = params["z"]
        self.e = params["e"]
        self.kb = params["kb"]
        self.T = params["T"]
        self.Na = params["Na"]
        self.Fa = params["Fa"]
        self.R = params["R"]
        self.D1 = params["D1"]
        self.D2 = params["D2"]
        self.Cb = params["Cb"]
        self.length = 1.0#params["length"]
        #d = 1.544E-6
        self.epsilon = params["epsilon"]
        self.coef =1 # 2 * self.e / ( self.kb * self.T )
        self.V0 =  self.coef * params["V0"]
        self.xi_0 = np.log(np.abs(np.tanh(self.V0/4)))
        self.g = np.tanh(self.xi_0/2)
        self.k = np.sqrt(2 * self.Cb * (self.z * self.Fa) ** 2 / (self.R * self.T * self.epsilon))
        self.d = 1/self.k
        self.xi_d = self.length
        
    def analytic_concentration(xi, params, r0 = 1E-5):

        z = params["z"]
        e = params["e"]
        kb = params["kb"]
        T = params["T"]
        Na = params["Na"]
        Fa = params["Fa"]
        R = params["R"]
        D1 = params["D1"]
        D2 = params["D2"]
        Cb = params["Cb"]
        #d = 1.544E-6
        epsilon = params["epsilon"]
        coef = 2*e/(kb*T)
        V0 =  coef * params["V0"]
        k = np.sqrt(2 * Cb * (z*Fa) ** 2 / (R * T * epsilon))

        def coth(x):
            return 1/np.tanh(x)

        def phi0(xi):
            return 2 * np.log(np.abs(np.tanh(0.5 * (xi-xi_0))))

        def C1(xi):
            return (1 / k) * np.exp( -phi0(xi) ) * ( xi - np.tanh(0.5 * (xi-xi_0)) - np.tanh( xi_0 / 2 ) )

        def C0(xi, s = 1):
            return Cb * np.exp( - s * phi0(xi) )

        def Cm0(xi):
            return Cb * np.exp( phi0(xi) )
        def C(xi):
            print("r value: " + str(r))
            return C0(xi) + r * C1(xi)

        def plot_concentration(xi, c, c0):
            plt.figure()
            plt.title('Analytic concentration to zero and first order \n in r', fontsize=16, fontweight='bold')
            plt.xlabel(r'$\xi = \kappa x$', fontsize=16)
            plt.ylabel('Molar Concentration ', fontsize=16)
            plt.plot(xi, c0, 'b', ms = 1, label='Zero order concentration')
            plt.plot(xi, c, 'g--', ms = 1, label='First order concentration')
            plt.legend(loc='upper right', shadow=True, fontsize='x-large').get_frame().set_facecolor('#00FFCC')
            plt.subplots_adjust(hspace=0.4)

            plt.grid(True, color= '#F2F2F2')
            plt.savefig('concentrations.eps', format='eps', dpi=1000, fontsize=16, fontweight='bold')
            plt.show()

        r = - r0

        xi_0 = np.log(np.abs(np.tanh(V0/4)))
        g = np.tanh(xi_0/2)


        c0 = C0(xi)
        c = C(xi)
        cm0 = Cm0(xi)

        #plot_concentration(xi, c, c0)
        file = open('../results/cp1r'+str(r0)+'.txt', 'w')
        file2 = open('../results/cp0r'+str(r0)+'.txt', 'w')
        file3 = open('../results/cm0r' + str(r0) + '.txt', 'w')

        for i in range(0, len(xi)):
            file.write(str(xi[i]/k*1E9) + "\t" +  str(c[i]) + "\n")
            file2.write(str(xi[i]/k*1E9) + "\t" +  str(c0[i]) + "\n")
            file3.write(str(xi[i]/k*1E9) + "\t" + str(cm0[i]) + "\n")

        file.close()
        return [xi/k*1E9, c, cm0]

    def analytic_potential(xi, params, r0 = 1e-5):
        def simpson(f, a, b, n = 100):
            h=(b-a)/n
            k=0.0
            x=a + h
            aux = int(n/2)
            for i in range(1,aux + 1):
                k += 4*f(x)
                x += 2*h

            x = a + 2*h
            for i in range(1,aux):
                k += 2*f(x)
                x += 2*h
            return np.array((h/3)*(f(a)+f(b)+k))

        def coth(x):
            return 1/np.tanh(x)

        def sech(x):
            return 1/np.cosh(x)



        def phi0(xi):
            return 2 * np.log(np.abs(np.tanh(0.5 * (xi-xi_0))))


        def G(xi):
            return np.log(np.abs(np.cosh(0.5 * (xi-xi_0))))

        def I_G(xi):
            I = simpson(G, xi, xi_d,)
            return I[0]

        def phiprime(xi):
            return xi ** 2 / 2 - 2 * g * xi + 2 * (2 * g - xi ) * coth(0.5*(xi-xi_0))

        def phiprime0(xi):
            return coth( 0.5 * (xi-xi_0) ) * ( sech( 0.5 * (xi-xi_0) ) ) ** 2

        def phi1(xi):
            aRes = []
            c = -phiprime(length)
            A = -2*xi_d/g + 2*g*xi_d - xi_d ** 2 / 2 - 2 * g 
            B = -(xi_d - 2/g)
            C = 3/2
            D = -2 * (xi_d - 2/g)
            E = 2
            response = -(A * (xi_d - xi) + B/2 * (xi_d ** 2 - xi **2) + C/3 * (xi_d**3 - xi **3) + D * np.log(np.abs(np.cosh((xi_d-xi_0)/2)/np.abs(np.cosh((xi-xi_0)/2)))) + E * (xi_d * np.log(np.abs(np.cosh((xi_d-xi_0)/2)))) - xi * np.log(np.abs(np.cosh((xi-xi_0)/2)))-I_G(xi))
            return response

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


        z = params["z"]
        e = params["e"]
        kb = params["kb"]
        T = params["T"]
        Na = params["Na"]
        Fa = params["Fa"]
        R = params["R"]
        coef = 2*e/(kb*T)
        V0 = coef * params["V0"]
        D1 = params["D1"]
        D2 = params["D2"]
        Cb = params["Cb"]
        length = params["length"]
        #d = 1.544E-6

        epsilon = params["epsilon"]
        #########################  PLOTTING SCRIPT  ###############################################
        k = np.sqrt(2 * Cb * (z*Fa) ** 2 / (R * T * epsilon))
        d = 1/k
        r = - r0

        xi_0 = np.log(np.abs(np.tanh(V0/4)))
        g = coth(xi_0/2)
        xi_d = xi[-1]

        p1 = 1/coef * phi(xi)
        p0 = 1/coef * phi0(xi)

        E1 = 1/coef * Efield(xi)
        E0 = 1/coef * Efield0(xi)

        #plot_potential(xi, p0, p1)

        file = open('../results/potential0.txt', 'w')
        file2 = open('../results/E-ana-r' + str(r0) + '.txt', 'w')
        file3 = open('../results/E1-ana-r' + str(r0) + '.txt', 'w')

        for i in range(0, len(xi)):
            file.write(str(xi[i]/k*1E9) + "\t" +  str(p0[i]) + "\n")
            file2.write(str(xi[i]/k*1E9) + "\t" + str(E0[i]) + "\n")
            file3.write(str(xi[i]/k*1E9) + "\t" + str(E1[i]) + "\n")

        file.close()

        file = open('../results/potential1-A'+str(r0)+'.txt', 'w')
        for i in range(0, len(xi)):
            file.write(str(xi[i]/k*1E9) + "\t" +  str(p1[i]) + "\n")

        file.close()
        file2.close()
        file3.close()

        return [xi/k*1E9, p1, E1]