# -*- coding: utf-8 -*-
"""
Created on Wed Jul 27 17:45:46 2016

@author: daniella
"""

def qdist(nsamples):
    
    import numpy as np
    import scipy
    #From Allen 2007
    gamma = 1.8 
    
    randomq = scipy.stats.powerlaw.rvs(gamma+1,size=nsamples)

#    plt.hist(test)
#    plt.hist(scipytest)
#    hy,hx = np.histogram(test)
#    y = max(hy)*x**1.8
#    plt.plot(x,y)

    return randomq
