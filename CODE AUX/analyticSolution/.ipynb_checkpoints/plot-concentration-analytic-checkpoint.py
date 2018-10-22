import matplotlib.pyplot as plt
import numpy as np
from system_definition import *

#########################  PLOTTING SCRIPT  ###############################################
#plt.style.use('seaborn')

N = 1E5
length = d()
step = np.float_(length/N)
x = np.arange(0, length, step)

## Concentration Points ##
cplus, cminus = conc_0order(phi0, x)
cplusI = conc_1order(phi0, x)

## PLOTS ##
fig1 = plt.figure(1, figsize=(16,12))

plt.subplot(2,1,1)
plt.title('Concentration close to the interfase', fontsize=20, fontweight='bold')
plt.ylabel('Molar Concentration ', fontsize=20)
plt.plot(x, cminus, 'b', ms = 1, label='Equilibrium Concentration of $SO_4^-$ ions')
plt.plot(x, cplus, 'r--', ms = 1, label='Equilibrium Concentration of $Cu^{+2}$ ions')
plt.legend(loc='upper right', shadow=True, fontsize='x-large').get_frame().set_facecolor('#00FFCC')

x_labels = [0, '$1 x 10^{-5}$','$2 x 10^{-5}$','$3 x 10^{-5}$','$4 x 10^{-5}$']
x_ticks = [0,0.000001 , 0.000002,0.000003 ,0.000004]
plt.xticks(x_ticks, x_labels)


plt.grid(True, color= '#F2F2F2')

plt.subplot(2,1,2)
#plt.title('Concentration close to the interfase', fontsize=16, fontweight='bold')
plt.xlabel('Distance to the interface (m)', fontsize=20)
plt.ylabel('Molar Concentration ', fontsize=20)
plt.plot(x, cplusI, 'g', ms = 1, label='Equilibrium Concentration $Cu^{+2}$ first order')
plt.plot(x, cplus, 'r--', ms = 1, label='Equilibrium Concentration $Cu^{+2}$ ion')
plt.legend(loc='upper right', shadow=True, fontsize='x-large').get_frame().set_facecolor('#00FFCC')
plt.subplots_adjust(hspace=0.4)

plt.grid(True, color= '#F2F2F2')

plt.xticks(x_ticks, x_labels)


## THIS COMMAND PRINTS THE FIGURE (COMPLETE) TO EPS
plt.savefig('concentrations.eps', format='eps', dpi=1000, fontsize=16, fontweight='bold')

plt.show()


