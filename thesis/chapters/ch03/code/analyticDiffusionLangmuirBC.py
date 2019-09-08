import numpy as np
from bisection import findRoot
import json

with open('model_parameters.json', 'r') as file:
    params = json.loads(file.read())


def C_an(xi, tau, LAMBDA_, params, tol = 1e-12):
    def cot(x):
        return np.cos(x)/np.sin(x)
    
    def csc(x):
        return 1/np.sin(x)
    
    def A(a):
        return 1 
    
    def P(a, xi):
        return  - 2 / (1 + kf * d / D ) * (a * np.tan(a/2) + a * ( csc(a) - 1 ) )/( a ** 2 * csc(a) ** 2 + a )*( np.cos( a * xi ) - cot(a) * np.sin( a * xi ) ) 
    
    total = np.zeros(len(xi))
    
    Cb = params['bulkConcentration']
    D = params['diffusionCoefficientCu']
    d = params['laminarFlowRegion']
    kf = params['reactionRate']
    i0 = params['reactionRate']
    ## Summing contributions
    for i in range(0,len(LAMBDA_)):
        a = LAMBDA_[i]
        total = total - ( 2 / (1 + kf * d / D) ) * np.exp( - ( a ** 2 ) * tau) *  (-kf * d / D + kf * d / D * a * csc(a) + a * np.tan(a/2) ) / ( kf * d / D + (a * csc(a)) ** 2) * ( np.cos( a * xi ) - cot(a) * np.sin( a * xi ) ) 
        
    C_ss = Cb * (1 + kf * d / D * xi )/( 1 + kf * d / D )
    print("Computation finished")
    return  C_ss + Cb * total



def findLambdas(params, N = 1000):

    def f(z, m): 
        # a = (2 * m + 1) * np.pi / 2 * z
        return np.tan((2 * m + 1) * np.pi / 2 * z) + D / (kf * d) * (2 * m + 1) * np.pi / 2 * z
        #return np.tan(z) + D / (kf * d) * z

    # note that tau = D t / xb ** 2, where t is the time
    Cb = params['bulkConcentration']
    D = params['diffusionCoefficientCu']
    d = params['laminarFlowRegion']
    kf = params['reactionRate']
    cond = True
    
    print("Starting iteration")
    er = 0.5e-15
    

    #print("Finding " + str(2*N) + " roots...")
    print("Finding " + str(N) + " roots...")
    LR = []
    LL = []
    #Forward
    for i in range(0,N):
        #Search for root of f (find lambda --a in this code)
        #a = (2 * i + 1) * np.pi / 2 + er
        #b = (2 * i + 3) * np.pi /2 - er
        a = 1 + er
        b = 3 - er
        z = findRoot(f, i, a, b)
        if (z == None):
            print("Found NoneType Lambda at N = "+str(i))
            break
        lam = (2 * i + 1 ) * np.pi / 2 * z
        LR.append(lam)
    LR = np.array(LR) 
    #backwards   
    #for i in range(0,N):
        #Search for root of f (find lambda --a in this code)
    #    a = -1 - er
    #    b = -3 + er
    #    z = findRoot(f, i, a, b)
    #    if (z == None):
    #        print("Found NoneType Lambda at N = "+str(i))
    #        break
    #    lam = (2 * i + 1 ) * np.pi / 2 * z
    #    LL.append(lam) 
    LL = -1 * np.array(LR)
    print("Done Computing lambdas")
    return LL, LR 