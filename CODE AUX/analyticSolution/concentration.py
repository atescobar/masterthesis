import matplotlib.pyplot as plt
import scipy.integrate as integrate
import numpy as np

def analytic_concentration(r0 = 1E-5):
    def coth(x):
        return 1/np.tanh(x)

    def phi0(xi):
        return 2 * np.log(np.abs(np.tanh(0.5 * (xi-xi_0))))

    def C1(xi):
        return (1 / k) * np.exp( -phi0(xi) ) * ( xi - np.tanh(0.5 * (xi-xi_0)) - np.tanh( xi_0 / 2 ) )

    def C0(xi):
        return C_b * np.exp( -phi0(xi) )

    def Cm0(xi):
        return C_b * np.exp( phi0(xi) )
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


    N = 1E4
    length = 20
    step = np.float_(length/N)
    xi = np.arange(0, length, step)
    C_b = 0.1
    k = 6.476E7
    
    r = - k * r0
    d = 1.544E-6
    V0 = - 77 * 0.15
    xi_0 = np.log(np.abs(np.tanh(V0/4)))
    g = np.tanh(xi_0/2)


    c0 = C0(xi)
    c = C(xi)
    cm0 = Cm0(xi)
    
    #plot_concentration(xi, c, c0)
    file = open('./results/cp1r'+str(r0)+'.txt', 'w')
    file2 = open('./results/cp0r'+str(r0)+'.txt', 'w')
    file3 = open('./results/cm0r' + str(r0) + '.txt', 'w')

    for i in range(0, len(xi)):
        file.write(str(xi[i]) + "\t" +  str(c[i]) + "\n")
        file2.write(str(xi[i]) + "\t" +  str(c0[i]) + "\n")
        file3.write(str(xi[i]) + "\t" + str(cm0[i]) + "\n")

    file.close()

#reactRate = [5E-3, 3E-3, 5E-3]
reactRate = [0, 8E-4, 7E-3]
for i in range(0,len(reactRate)):
    analytic_concentration(reactRate[i])

