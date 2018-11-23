import numpy as np
import matplotlib.pyplot as plt


def C_an(xi, tau):
    # note that tau = D t / xb ** 2, where t is the time
    Cb = 100
    total = 0
    tolMax  = 1e-3* np.ones(len(xi))
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
            contrib = np.abs(actual)/np.abs(previo)
            for i in range(0, len(xi)):
                if ( contrib[i] > tolMin[i]):
                    cond = True
                if ( contrib[i] < tolMax[i]):
                    cond = False
                        
        previo = actual
    print("Series truncation at n = " + str(m) + " with average local error " + str(contrib.sum()/len(contrib) * 100) + "%")
    return  Cb - ( (4 * Cb / np.pi ) * total ) 







