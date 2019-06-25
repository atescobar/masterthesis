import numpy as np
import json

with open('model_parameters.json', 'r') as file:
    params = json.loads(file.read())
    
def C_an(xi, tau, r):
    # note that tau = D t / xb ** 2, where t is the time
    Cb = params['bulkConcentration']
    D = params['diffusionCoefficientCu']#Diffusion Coefficient
    d = 1e-9 #Laminar flow sheet
    totalA = 0
    totalB = 0
    tolMax  = 1e-10 
    tolMin = tolMax*1e-2
    m = 0
    cond = True
    previo = 10 * np.ones(len(xi))
    
    while(cond):
        Am = ((-1) ** m ) / ( 2 * m + 1 ) * np.exp( - ( ( 2 * m + 1 ) * ( np.pi / 2 ) ) ** 2 * tau ) * np.cos( (2 * m + 1) * np.pi / 2 * xi )
        Bm = 1 / ( 2 * m + 1 ) ** 2 * np.exp( - ( ( 2 * m + 1 ) * ( np.pi / 2 ) ) ** 2 * tau ) * np.cos( (2 * m + 1) * np.pi / 2 * xi )
        actualA = Am
        actualB = Bm
        totalA += actualA
        totalB += actualB
        m = m + 1
        
        if (m > 1):
            contrib = (actualA)/previo
            for i in range(0, len(xi)-1):
#                if ( np.abs(contrib[i]) > tolMin):
#                    cond = True
                if ( np.abs(max(contrib)) < tolMax):
                    
                    cond = False
                        
        previoA = actualA
    print("Series truncation at n = " + str(m) + " with maximnum local error " + str(np.amax(np.abs(contrib))) + "%")
    return  Cb - ( (4 * Cb / np.pi ) * totalA + ( r * d / D ) * (np.ones(len(xi)) - xi - ( 8 / np.pi ** 2 ) * totalB )) 







