import numpy as np


def C_an(xi, tau, params):
    # note that tau = D t / xb ** 2, where t is the time
    Cb = params['bulkConcentration']
    total = 0
    tolMax  = 1e-10 
    tolMin = tolMax*1e-2
    m = 0
    cond = True
    previo = 10 * np.ones(len(xi))
    
    while(cond):
        Am = ((-1) ** m ) / ( 2 * m + 1 ) * np.exp( - ( ( 2 * m + 1 ) * ( np.pi / 2 ) ) ** 2 * tau ) 
        actual = Am * np.cos( (2 * m + 1) * np.pi / 2 * xi )
        total += actual
        m = m + 1
        
        if (m > 1):
            contrib = (actual)/previo
            for i in range(0, len(xi)-1):
#                if ( np.abs(contrib[i]) > tolMin):
#                    cond = True
                if ( np.abs(max(contrib)) < tolMax):
                    
                    cond = False
                        
        previo = actual
    print("Series truncation at n = " + str(m) + " with maximnum local error " + str(np.amax(np.abs(contrib))) + "%")
    return  Cb - ( (4 * Cb / np.pi ) * total ) 







