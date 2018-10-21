
import numpy
import matplotlib.pyplot as plt
#from analyticSolution.concentration import analytic_concentration as aConc
#from analyticSolution.potential import analytic_potential as aPot
#from numericalSolution.diffusionReactionSystem import diffusionReactionSystem  as dRS

def plotOne(x, F, title, xLabel= "", yLabel = "", FLabel = ""):
    plt.figure(1)
    plt.title(title, fontsize=16, fontweight='bold')
    plt.xlabel(r'$\xi = \kappa x$', fontsize=16)
    plt.ylabel(r'Dimentionless potential ', fontsize=16)
    plt.plot(xi, phi, 'g', label = "Potential")
    plt.axis([xStart, xStop, -20, -0.0])
    plt.legend()
    ## THIS COMMAND PRINTS THE FIGURE (COMPLETE) TO EPS
    plt.savefig(shortname + '.eps', format='eps', dpi=1000, fontsize=16, fontweight='bold')
    plt.show()

def plotTwo(xf, F, xg, G, shortname = "1", title = "", xLabel= "", yLabel = "", Flabel = "", Glabel = ""):
    plt.figure(1)
    plt.title(title, fontsize=16, fontweight='bold')
    plt.xlabel(xLabel, fontsize=16)
    plt.ylabel(yLabel, fontsize=16)
    plt.plot(xf, F, 'g', label = Flabel)
    plt.plot(xg, G, 'g', label = Glabel)
    #plt.axis([xStart, xStop, -20, -0.0])
    plt.legend()
    ## THIS COMMAND PRINTS THE FIGURE (COMPLETE) TO EPS
    plt.savefig(shortname + '.eps', format = 'eps', dpi = 1000, fontsize=16, fontweight='bold')
    plt.show()

def getResults(file):
    results = open(file, 'r')
    xi = []
    f = []
    for line in results:
        u, v = line.split()
        xi.append(u)
        f.append(v)
    return xi, f


#reactRate = [5E-3, 3E-3, 5E-3]
reactRate = [5E-8, 3E-8, 1E-8]
r0 = reactRate[0]
r1 = reactRate[1]
r2 = reactRate[2]


filePhiN = './results/potential-num-r'+str(r0)+'.txt'
filePhiN1 = './results/potential-num-r'+str(r1)+'.txt'
filePhiN2 = './results/potential-num-r'+str(r2)+'.txt'

fileCpN = './results/cp-num-r'+str(r0)+'.txt'
fileCmN = './results/cm-num-r'+str(r0)+'.txt'
filePhiA = './results/potential1-A'+str(r0)+'.txt'
fileCmA = './results/cm0r'+str(r0)+'.txt'
fileCpA = './results/cp1r'+str(r0)+'.txt'


fileEN = './results/E-num-r' + str(r0) + '.txt'
fileEN1 = './results/E-num-r' + str(r1) + '.txt'
fileEN2 = './results/E-num-r' + str(r2) + '.txt'

fileEA = './results/E-ana-r' + str(r0) + '.txt'
fileEA1 = './results/E-ana-r' + str(r1) + '.txt'
fileEA2 = './results/E-ana-r' + str(r2) + '.txt'


xiNpn, phiN = getResults(filePhiN)
xiNpn1, phiN1 = getResults(filePhiN1)
xiNpn2, phiN2 = getResults(filePhiN2)

xiNE, EN = getResults(fileEN)
xiNE1, EN1 = getResults(fileEN1)
xiNE2, EN2 = getResults(fileEN2)

xiAE, EA = getResults(fileEA)
xiAE1, EA1 = getResults(fileEA1)
xiAE2, EA2 = getResults(fileEA2)


xiNcpn, CpN = getResults(fileCpN)
xiNcmn, CmN = getResults(fileCmN)
xiAcpa, CpA = getResults(fileCpA)
xiAcma, CmA = getResults(fileCmA)
xiApa, phiA = getResults(filePhiA)

##  PLOT NUMERIC ONLY
#plotTwo(xiNcpn, CpN, xiNcmn, CmN, 'conN' + append, 'Numeric Concentration', r'$\kappa \delta$', r'Molar Concentration', r'C_+', r'C_-')
#plotOne(xiNpn, phiN, 'phiN'+ append, 'Numeric Potential', r'$\kappa \delta$', r'Dimensionless Potentail', r'\Phi')


##  PLOT ANALYTIC ONLY
#plotTwo(xiAcpa, CpA, xiAcmn, CmA, 'conA' + append, 'Analytic Concentration', r'$\kappa \delta$', r'Molar Concentration', r'C_+', r'C_-')
#plotOne(xiApa, phiA, 'phiA'+ append, 'Analytic Potential', r'$\kappa \delta$', r'Dimensionless Potentail', r'\Phi')

##  PLOT COMPARISON ANALYTIC-NUMERIC
#plotTwo(xiApa, phiA, xiNpn, phiN, 'potAN'+ append, 'Comparison Of Analytic and Numeric Potential', r'$\kappa \delta$', r'Dimentionless Potential', r'\Phi_A', r'\Phi_N')
#plotTwo(xiAcpa, CpA, xiNcpn, CpN, 'CpAN'+ append, 'Comparison Of Analytic and Numeric Concentration', r'$\kappa \delta$', r'Molar Concentration', r'C_+', r'C_-')

## Plot Comparison Conc
plt.figure(1)
plt.title("Concentration comparison", fontsize=16, fontweight='bold')
plt.xlabel("Adimentional distance", fontsize=16)
plt.ylabel("Molar Concentration", fontsize=16)
plt.plot(xiNcpn, CpN, 'r--', label = r'Numeric $C_+$')
plt.plot(xiAcpa, CpA, 'r', label = r'Analytic $C_+$')
plt.plot(xiNcmn, CmN, 'b--,', label = r'Numeric $C_{-}$')
plt.plot(xiAcma, CmA, 'b', label = r'Analytic $C_{-}$')
plt.ylim(0.0 ,1.0)
#plt.axis([xStart, xStop, -20, -0.0])
plt.ylim(ymax=1.0)
plt.legend()
## THIS COMMAND PRINTS THE FIGURE (COMPLETE) TO EPS
plt.savefig("concentration-comparison" + '.eps', format = 'eps', dpi = 1000, fontsize=16, fontweight='bold')
plt.show()

# Plot Comparison Potential
plt.figure(2)
plt.title("Comparison Potential", fontsize=16, fontweight='bold')
plt.xlabel("Adimentional Distance", fontsize=16)
plt.ylabel("Dimentionless Potential", fontsize=16)
#plt.plot(xiNpn1, phiN1, 'g', markersize=1,  label = "Numeric Potential r = " + str(r1))
plt.plot(xiNpn, phiN, 'r--', markersize=1,  label = "Numeric Potential r = " + str(r0))
#plt.plot(xiNpn2, phiN2, 'bx', markersize=1, label = "Numeric Potential r = " + str(r2))
#plt.plot(xiApa, phiA, 'b', label = "Analytic Potential")
#plt.axis([xStart, xStop, -20, -0.0])
plt.legend()
## THIS COMMAND PRINTS THE FIGURE (COMPLETE) TO EPS
plt.savefig("potential-comparison" + '.eps', format = 'eps', dpi = 1000, fontsize=16, fontweight='bold')
plt.show()

plt.figure(3)
plt.title("Numeric Potential", fontsize=16, fontweight='bold')
plt.xlabel("Adimentional Distance", fontsize=16)
plt.ylabel("Dimentionless Potential", fontsize=16)

plt.plot(xiNpn, phiN, 'g', markersize=2, label = "Numeric Potential r = " + str(r0))
#plt.plot(xiNpn1, phiN1, 'b--', markersize=1, label = "Numeric Potential r = " + str(r1))
#plt.plot(xiNpn2, phiN2, 'rx', markersize=1,label = "Numeric Potential r = " + str(r2))

#plt.axis([xStart, xStop, -20, -0.0])
plt.legend()
## THIS COMMAND PRINTS THE FIGURE (COMPLETE) TO EPS
plt.savefig("potential-numeric-r = " + str(r0) + '.eps', format = 'eps', dpi = 1000, fontsize=16, fontweight='bold')
plt.show()



plt.figure(4)
plt.title("Numeric Electric Field", fontsize=16, fontweight='bold')
plt.xlabel("Adimentional Distance", fontsize=16)
plt.ylabel("Dimentionless Electric Field", fontsize=16)
plt.plot(xiNE, EN, 'r', markersize=1,  label = "Numeric Electric Field r = " + str(r2))
#plt.plot(xiNE1, EN1, 'g--', markersize=1, label = "Numeric Electric Field r = " + str(r0))
#plt.plot(xiNE2, EN2, 'bx', markersize=1, label = "Numeric Electric Field r = " + str(r1))

#plt.axis([xStart, xStop, -20, -0.0])
plt.legend()
## THIS COMMAND PRINTS THE FIGURE (COMPLETE) TO EPS
plt.savefig("Efield-numeric-r = " + str(r0) + '.eps', format = 'eps', dpi = 1000, fontsize=16, fontweight='bold')
plt.show()


plt.figure(5)
plt.title("Electric Potential Results", fontsize=16, fontweight='bold')
plt.xlabel("Adimentional Distance", fontsize=16)
plt.ylabel("Dimentionless Potential", fontsize=16)
plt.plot(xiApa, phiA, 'g', markersize=1,label = "Analytic Potential r = " + str(r0))
plt.plot(xiNpn, phiN, 'b--', markersize=1, label = "Numeric Potential r = " + str(r0))

#plt.axis([xStart, xStop, -20, -0.0])
plt.legend()
## THIS COMMAND PRINTS THE FIGURE (COMPLETE) TO EPS
plt.savefig("potential-results-r = " + str(r0) + '.eps', format = 'eps', dpi = 1000, fontsize=16, fontweight='bold')
plt.show()



plt.figure(6)
plt.title("Electric Field Results", fontsize=16, fontweight='bold')
plt.xlabel("Adimentional Distance", fontsize=16)
plt.ylabel("Dimentionless Electric Field", fontsize=16)
plt.plot(xiNE, EN, 'b--', markersize=2,  label = "Numeric E(x)r = " + str(r0))
plt.plot(xiAE, EA, 'g', markersize=2, label = "Analytic E(x) r = " + str(r0))
plt.ylim(0, 10)
#plt.axis([xStart, xStop, -20, -0.0])
plt.legend()
## THIS COMMAND PRINTS THE FIGURE (COMPLETE) TO EPS
plt.savefig("Efield-results-r = " + str(r0) + '.eps', format = 'eps', dpi = 1000, fontsize=16, fontweight='bold')
#plt.show()


