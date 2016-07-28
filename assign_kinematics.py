import numpy as np
from scipy import interpolate
from astropy.table import Table
import matplotlib.pyplot as plt

def uvw(tau0):

    # Taken from Aumer & Binney 2009
    # Assigns the dispersions for each velocity component
    
    # First do the U velocity component
    beta, tau1, v10 = 0.307, 0.001, 41.899
    sigU = v10 * ( (tau0 + tau1) / (10 + tau1) ) ** beta

    # Now do the V velocity component
    beta, tau1, v10 = 0.430, 0.715, 28.823
    sigV = v10 * ( (tau0 + tau1) / (10 + tau1) ) ** beta

    # Now do the W velocity component
    beta, tau1, v10 = 0.445, 0.001, 23.831
    sigW = v10 * ( (tau0 + tau1) / (10 + tau1) ) ** beta
    
    return sigU, sigV, sigW
    
def uvw_ave(tau0):

    'Assign average velocity components'
    'based on age'

    u = 0
    w = 0

    # Fit to Gontcharov 2012
    x, y = [0,5], [0,-19]
    z = np.polyfit(x, y, 1)
    p = np.poly1d(z)
    v = p(tau0)

    return u,v,w

def uvw_ave2(sigU):

    'Assign average velocity components'
    'based on velocity dispersion'

    u = 0
    w = 0

    # Stromberg asymmetric drift equation
    v = -1 * sigU**2 / 74.

    return u,v,w


def test():

    u,v,w = uvw(0.5)
    #um,vm,wm = uvw_ave(10)
    #print um,vm,wm
    um,vm,wm = uvw_ave2(u)
    #print um,vm,wm
    Us = np.random.normal(um, u, 10000)
    Vs = np.random.normal(vm, v, 10000)
    Ws = np.random.normal(wm, w, 10000)

    u1,v1,w1 = uvw(3)    
    #um1,vm1,wm1 = uvw_ave(1)
    #print um1,vm1,wm1 
    um1,vm1,wm1 = uvw_ave2(u1)
    #print um1,vm1,wm1
    Us1 = np.random.normal(um1, u1, 10000)
    Vs1 = np.random.normal(vm1, v1, 10000)
    Ws1 = np.random.normal(wm1, w1, 10000)

    u2,v2,w2 = uvw(10)
    #um,vm,wm = uvw_ave(10)
    #print um,vm,wm
    um2,vm2,wm2 = uvw_ave2(u2)
    #print um,vm,wm
    Us2 = np.random.normal(um2, u2, 10000)
    Vs2 = np.random.normal(vm2, v2, 10000)
    Ws2 = np.random.normal(wm2, w2, 10000)

    f, axarr = plt.subplots(2, sharex=True, figsize=(6,8) )

    axarr[0].scatter(0, 0, s = 10, c='r', edgecolors='None', label='500 Myr', zorder=-10)
    axarr[0].scatter(0, 0, s = 10, c='b', edgecolors='None', label='3 Gyr', zorder=-11)
    axarr[0].scatter(0, 0, s = 10, c='k', edgecolors='None', label='10 Gyr', zorder=-12)
    axarr[0].legend(loc=2, scatterpoints=1)
    axarr[0].scatter(Vs, Us, alpha=0.5, s=1, c='r', edgecolors='None', zorder=10)
    axarr[0].scatter(Vs1, Us1, alpha=0.5, s=1, c='b', edgecolors='None', zorder=9)
    axarr[0].scatter(Vs2, Us2, alpha=0.5, s=1, c='k', edgecolors='None', zorder=8)

    
    axarr[1].scatter(Vs, Ws, alpha=0.5, s=1, c='r', edgecolors='None', zorder=10)
    axarr[1].scatter(Vs1, Ws1, alpha=0.5, s=1, c='b', edgecolors='None', zorder=9)
    axarr[1].scatter(Vs2, Ws2, alpha=0.5, s=1, c='k', edgecolors='None', zorder=8)
    axarr[0].axhline(0, ls='--', c='k', alpha=0.5)
    axarr[0].axvline(0, ls='--', c='k', alpha=0.5)
    axarr[1].axhline(0, ls='--', c='k', alpha=0.5)
    axarr[1].axvline(0, ls='--', c='k', alpha=0.5)
    axarr[0].set_ylabel('U (km/s)')
    axarr[1].set_xlabel('V (km/s)')
    axarr[1].set_ylabel('W (km/s)')
    axarr[0].set_aspect('equal', 'datalim')
    axarr[1].set_aspect('equal', 'datalim')
    plt.tight_layout()
    plt.show()





