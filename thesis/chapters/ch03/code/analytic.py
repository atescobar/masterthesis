import numpy as np
import matplotlib.pyplot as plt


def C_an(xi, tau):
    # note that tau = D t / xb ** 2, where t is the time
    Cb = 100
    total = 0
    n = 100000
    
    for m in range(0, n):
        Am = ((-1) ** m ) / ( 2 * m + 1 ) * np.exp( - ( ( 2 * m + 1 ) * ( np.pi / 2 ) ) ** 2 * tau ) 
        total += Am * np.cos( (2 * m + 1) * np.pi / 2 * xi )
        
    return  Cb - ( (4 * Cb / np.pi ) * total ) 







