import numpy as np


def C_an(xi, tau, Cb):
    # note that tau = D t / xb ** 2, where t is the time
    F = 96485.3415
    kf = 1e1
    i0 = 150
    i1 = 0 #i0/(F * kf)
    total = np.zeros(len(xi))
    tolMax  = 1e-2
    tolMin = tolMax*1e-2
    m = 1
    cond = True
    previo = 10 * np.ones(len(xi))
    while(cond):
        Am = (((-1) ** m ) * Cb - i1) / m  * np.exp( - ( m * np.pi ) ** 2 * tau ) 
        actual = Am * np.sin( m * np.pi * xi )
        total += actual
    
        if (m > 2):
            contrib = Am/Am_previo
            for i in range(0, len(xi)-1):
                if ( np.abs(contrib) < tolMax):
                    cond = False
        
        m = m + 1
        Am_previo = Am
        previo = actual
    print("Series truncation at n = " + str(m) + " with maximnum local error " + str(np.abs(contrib)) + "%")
    return  i1 + (Cb - i1) * xi + (2 / np.pi ) * total







