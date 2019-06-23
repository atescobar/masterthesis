class Model():
    def __init__(self, params):
        #Model Parameters
        self.Cb = params['bulkConcentration']
        self.D = params['diffusionCoefficientScale']
        self.d = params['laminarFlowRegion']
        self.kf = params['reactionRate']
        self.z = params['z']
        self.F = params['Fa']
        self.R = params['R']
        self.T = params['T']
        self.epsilon = params['epsilon']
        self.kappa =  np.sqrt(( ( self.z * self.F  ) ** 2 * self.Cb ) / ( self.epsilon * self.R * self.T ) )
        self.Psi0 = self.z * self.F * params['V0'] / ( self.R * self.T )
        self.D1 = self.D * params["D1"]
        self.D2 = self.D * params["D2"]
        self.N = 100000
        self.M = 200
        self.xi = np.linspace(0,params["length"], self.M)
        self.tau = np.linspace(0,1, self.N) #shape is N+1

        #Grid Parameters
        self.dtau = 1/(self.N)  # N Partitions
        self.dxi = 1/(self.M) # N Partitions 
        self.a1 = self.dtau / self.dxi ** 2 
        self.a2 = self.dtau / self.dxi ** 2 * self.D2/self.D1

        #Plotting parameters
        self.imageName = 'complete-diffusion-nernst'
        
    def build(self):
        M = self.M
        N = self.N
        a1 = self.a1
        a2 = self.a2
        Psi0 = self.Psi0
        kappa = self.kappa
        kf = self.kf
        dxi = self.dxi
        D1 = self.D1 
        D2 = self.D2
        # Define the coefficient matrix
        g1 = 1 / ( 1 + kf * dxi / ( D1 * kappa ) +  Psi0)
        di1 = ( 1 - 2 * a1 ) * np.ones(M-2)
        di1[0] = ( 1 - 2 * a1 + a1 * g1 )
        A1 = diags(np.array([ a1 * np.ones(M-3), di1, a1 * np.ones(M-3)]), [-1, 0, 1], shape=(M-2, M-2)).toarray()

        g2 = 1 / ( 1 - Psi0 )
        di2 = ( 1 - 2 * a2 ) * np.ones(M-2)
        di2[0] = ( 1 - 2 * a2 + a2 * g2 )
        A2 = diags(np.array([ a2 * np.ones(M-3), di2, a2 * np.ones(M-3)]), [-1, 0, 1], shape=(M-2, M-2)).toarray()

        B1 = np.zeros([M-2, M-2])
        B2 = np.zeros([M-2, M-2])

        D0 = diags(np.array([ np.ones(M-3), -2 * np.ones(M-2), np.ones(M-3)]), [-1, 0, 1], shape=(M-2, M-2)).toarray()
        Dinv = np.asarray(np.linalg.inv(D0))

        b1 = np.zeros(M-2)
        b1[-1] = a1 
        b2 = np.zeros(M-2)
        b2[-1] = a2

        bPsi = np.zeros(M-2)
        bPsi[0] = Psi0

        def B(s, Psi, n):
            diag =  (Psi[n, 1:M-1 ] - Psi[n, 0:M-2 ])
            diag2 =  (Psi[n, 1:M-2 ] - Psi[n, 2:M-1 ])
            PsiMatrix = diags(np.array([ diag , diag2 ]), [0, 1], shape=(M-2, M-2)).toarray()
            if s == 1:
                return  -1 * a1 * PsiMatrix
            if s == -1:
                return a2 * PsiMatrix
        # Set up initial conditions for C

        rho1 = np.zeros([N, M])
        rho2 = np.zeros([N, M])
        Psi = np.zeros([N, M])

        rho1[0, :] = 0
        rho1[0, -1] = 1    

        rho2[0, :] = 0
        rho2[0, -1] = 1  

        Psi[0, :] = 0
        Psi[0, 0] = Psi0

        #Starting iteration
        for n in range(0, N-1):

             # Update border condition
            g1 = 1 / ( 1 + kf * dxi / ( D1 * kappa ) - (Psi[n,1]- Psi0))
            A1[0,0] = ( 1 - 2 * a1 + g1 * a1 )

            g2 = 1 / ( 1 + (Psi[n, 1] - Psi[n,0]))
            A2[0,0] = ( 1 - 2 * a2 + g2 * a2 )


            rho1[n+1, 1:M-1] = np.matmul(A1, rho1[n, 1:M-1])  + b1 + np.matmul(B(1, Psi, n), rho1[n, 1:M-1])
            rho1[n+1, 0] = g1 * rho1[n+1, 1]
            rho1[n+1, -1] = 1

            rho2[n+1, 1:M-1] = np.matmul(A2, rho2[n, 1:M-1]) + b2 + np.matmul(B(-1, Psi, n), rho2[n, 1:M-1]) 
            rho2[n+1, 0] = g2 * rho2[n+1, 1]
            rho2[n+1, -1] = 1

            Psi[n+1, 1:M-1] = np.matmul(Dinv, dxi * (rho2[n+1, 1:M-1] - rho1[n+1, 1:M-1]) -bPsi )
            Psi[n+1, 0] = Psi0
            Psi[n+1, -1] = 0
            
        print("Build Complete")
        self.rho1 = rho1
        self.rho2 = rho2
        self.Psi = Psi
        
    def remove_points(self, A, n):
        #n is the number of steps to skip
        if n >= 4:
            A = np.delete(A, [1, 2, 3])

        for i in range(0,int(len(A)/4)):
            index = i+n
            A = np.delete(A, [index-2, index-1, index])
        return A
        #Cm is the imported analytical solution


    def plot(self, t, imageName='complete-diffusion-nernst'):
        self.imageName = imageName
        Cb = self.Cb
        dtau = self.dtau
        C1 = Cb * self.rho1
        C2 = Cb * self.rho2
        phi = self.R * self.T * self.Psi / (self.z * self.F)
        kappa = self.kappa
        D1 = self.D1 
        mw = 4
        fs = 24
        skip = 4
        nanometerScale = abs(int(math.log10(1/self.kappa)))/9


        #xi2 = remove_points(self.xi, skip) # this is done to avoid cluttering of numeric points over the analytic solution
        xi2 = self.xi # this is done to avoid cluttering of numeric points over the analytic solution
        xi2 = xi2 * nanometerScale #change the scale of the scale to nanometer

        fig, ax1 = plt.subplots(figsize=(20,16))


        color = 'tab:red'
        ax1.tick_params(axis='y', labelcolor=color)
        ax2 = ax1.twinx() 
        color = 'tab:blue'

        plt.title('Comparing Numeric Solution To The Diffusion Reaction And Analytic Solution \n To The Diffusion-Only Problem', fontsize=fs, fontweight='bold')

        n = int(t/dtau)
        #print("Time: " + str(t / (D1 * kappa ** 2)) + "s")
        #ax1.plot(xi2, remove_points(C1[n, :], skip), 'g^', markersize=mw, label=r'$C_+$,  $\tau ='+str(t)+'$')
        #ax1.plot(xi2, remove_points(C2[n, :], skip), 'r^', markersize=mw, label=r'$C_-$,  $\tau ='+str(t)+'$')
        #ax1.legend(loc='upper left', fontsize = fs-4)
        #ax2.plot(xi2, remove_points(phi[n, :], skip), 'b^', markersize=mw, color='tab:blue', label=r'$\phi$,  $\tau ='+str(t)+'$')
        
        ax1.plot(xi2, C1[n, :], 'g^', markersize=mw, label=r'$C_+$,  $\tau ='+str(t)+'$')
        ax1.plot(xi2, C2[n, :],'r^', markersize=mw, label=r'$C_-$,  $\tau ='+str(t)+'$')
        ax1.legend(loc='upper left', fontsize = fs-4)
        ax2.plot(xi2, phi[n, :], 'b^', markersize=mw, color='tab:blue', label=r'$\phi$,  $\tau ='+str(t)+'$')

        ax1.set_xlabel(r'Distance from the interface plate (nm)', fontsize=fs)
        ax1.set_ylabel(r'Molar Concentration', fontsize=fs)
        ax2.set_ylabel(r'Electric Potential (V)', fontsize=fs)
        ax2.tick_params(axis='y', labelcolor=color)

        #plt.text(0.9, 40, r'Reaction Rate', fontsize = 14, color = 'black')
        #plt.text(0.9, 37, r'$r = - 1.5 \times 10 \frac{A}{m^2}$', fontsize = 14, color = 'black')
        plt.legend(loc='upper right', fontsize = fs-4)
        if(len(imageName) > 0):
            plt.savefig('../../../img/'+ self.imageName +'.eps', format='eps', dpi=1000, fontsize=16, fontweight='bold')
            
            
        fig.tight_layout()  # otherwise the right y-label is slightly clipped

    
        plt.show()



            
    

