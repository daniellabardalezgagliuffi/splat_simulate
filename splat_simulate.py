"""
Authors:
    Christian Aganze
    Daniella Bardalez Gagliuffi
    Caleb Choban
    Chris Theissen
"""

import splat
import numpy as np
import os
from astropy.table import Table
import scipy.interpolate as interpolate
import matplotlib.pyplot as plt
from scipy.special import erfc
from scipy.stats import powerlaw
import pandas as pd
import warnings

def test():

    #f, axarr = plt.subplots(4, figsize=(6,8) )

    #density = densities()
    #ages = ages(10000, [1,10])
    #masses = generate_masses([.01,.1])
    #starTable = make_star(masses,ages)

    u,v,w = uvw(0.5)
    um,vm,wm = uvw_ave2(u)
    Us = np.random.normal(um, u, 10000)
    Vs = np.random.normal(vm, v, 10000)
    Ws = np.random.normal(wm, w, 10000)

    u1,v1,w1 = uvw(3)    
    um1,vm1,wm1 = uvw_ave2(u1)
    Us1 = np.random.normal(um1, u1, 10000)
    Vs1 = np.random.normal(vm1, v1, 10000)
    Ws1 = np.random.normal(wm1, w1, 10000)

    u2,v2,w2 = uvw(10)
    um2,vm2,wm2 = uvw_ave2(u2)
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

def make_star(mass,age,model='burrow'):
    if np.min(mass) < 0.0005:
        raise NameError('Mass below minimum mass of 0.0005Msun')
    if np.max(mass) > 0.1 and model=='baraffe':
       warnings.warn('Mass above maximum mass of 0.1Msun for Baraffe 2003. Using Burrows 1997 instead.')  
       model = 'burrows'
    if np.min(mass) > 0.2:
        raise NameError('Mass above maximum mass of 0.2Msun for Burrows 1997')  
        
    model = splat.bdevopar.ReadModel(model)
    star_list = splat.bdevopar.Parameters(model, masses=mass, ages=age)
    return star_list

# This is broken Daniella.
def bad_make_star(mass, age, model='Burrows97'):
    
    '''
    Calculates stellar properties such as Teff, radii, logg and logL using evolutionary models
    '''

    if np.min(mass) < 0.0005:
        raise NameError('Mass below minimum mass of 0.0005Msun')

    if np.max(mass) > 0.1 and model=='Baraffe03':
       warnings.warn('Mass above maximum mass of 0.1Msun for Baraffe 2003. Using Burrows 1997 instead.')  
       model = 'Burrows97'

    if np.min(mass) > 0.2:
        raise NameError('Mass above maximum mass of 0.2Msun for Burrows 1997')    

    if model == 'Burrows97':
        #0.0005 - 0.2 Msun
        burrows = pd.read_pickle("burrows97.pickle")
        allages = burrows["Age (Gyr)"]
        allmasses = burrows["M/Ms"]
        teff = burrows["Teff"]
        radius = burrows["R/Rs"]
        logg = burrows["logg(cgs)"]
        logL = burrows["logL/Ls"]
        
    if model == 'Baraffe03':
        #0.0005 - 0.1 Msun
        baraffe = pd.read_pickle("baraffe03.pickle")
        allages = baraffe["Age (Gyr)"]
        allmasses = baraffe["M/Ms"]
        teff = baraffe["Teff"]
        radius = baraffe["R/Rs"]
        logg = baraffe["logg(cgs)"]
        logL = baraffe["logL/Ls"]

    interpteff = interpolate.interp2d(allages,allmasses,teff,kind='linear')
    interprad = interpolate.interp2d(allages,allmasses,radius,kind='linear')
    interplogg = interpolate.interp2d(allages,allmasses,logg,kind='linear')
    interplogL = interpolate.interp2d(allages,allmasses,logL,kind='linear')  
        
    mass = np.array(mass).flatten()
    age = np.array(age).flatten()
    
    newteff = np.array([interpteff(i,j) for i,j in zip(age,mass)])
    newrad = np.array([interprad(i,j) for i,j in zip(age,mass)])
    newlogg = np.array([interplogg(i,j) for i,j in zip(age,mass)])
    newlogL = np.array([interplogL(i,j) for i,j in zip(age,mass)])
    
    stardict = {'Teff (K)':newteff,'Radius (Rs)':newrad, 'log g':newlogg, 'log L':newlogL}    
    
    starTable = Table(stardict)
    
    return starTable

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


def ages(numstars, agerange=[0,10]):
    low, high = agerange
    ages = np.random.uniform(low, high, size=numstars)
    return ages

def densities(lffile = None, colors = None, interp = True, size = 10000):

    '''
    Reads in a file with a luminosity function and returns stellar density
    '''

    if lffile is None:
        lffile = os.path.dirname(os.path.realpath("__file__")) + \
                 '/LFs/Reyle_2010_J_LF.csv'

    # Check if the file exists
    if os.path.isfile(lffile) is False:
        raise LookupError("File %s does not exist"%lffile)

    # Read in the luminosity function
    t       = Table.read(lffile, format='ascii.csv')

    # Pull the values from the table
    AbsMag  = t['M_']
    PhiMean = t['Phi_Mean']
    PhiMin  = PhiMean - t['Phi_Minus']
    PhiMax  = PhiMean + t['Phi_Plus']

    # Interpolate the LF
    if interp is True: # Currently the only version that is implemented
        phiave  = interpolate.InterpolatedUnivariateSpline( AbsMag, PhiMean, k=1)
        philow  = interpolate.InterpolatedUnivariateSpline( AbsMag, PhiMin,  k=1)
        phihigh = interpolate.InterpolatedUnivariateSpline( AbsMag, PhiMax,  k=1)

    # Get the absolute magnitude range
    M1, M2 =  min(AbsMag), max(AbsMag)

    # Integrate the density function
    # First we pull a random phi from a triangular function.
    # 'size' defines how fine you want the interpolation grid to be 
    # (bigger is better, but computationally expensive)
    phi       = np.array( [ np.random.triangular( philow(x), phiave(x), phihigh(x), size ) \
                            for x in np.linspace(M1, M2, size)  ] )
    densities = np.array( [ np.trapz( x = np.linspace(M1, M2, size), y = phi[:,i] ) * 1e-3 \
                            for i in range(size) ] )

    return densities
    
# Chabrier CDF taken from Jumper & Fisher (2013) w/ params from Chabrier (2005)
def _cdf_lowmass_Chabrier(mass):
    A_s = 0.724;  sig = 0.55; m_c = 0.2;
    return A_s * sig * np.sqrt(np.pi/2.) * erfc((np.log10(m_c) - np.log10(mass)) / (sig * np.sqrt(2.0)))
def _cdf_highmass_Chabrier(mass):
    B_s = 0.323; x = 1.35; m_o = 1.; C = 0.896;
    return C + (B_s / np.log(10.)) * (np.power(m_o,-x) / x) * (1. - np.power((mass / m_o),-x))

# Using Kroupa (2001) w/ CDF calculated by hand
# k values used to connect each part of the Kroupa IMF
def _cdf_lowmass_Kroupa(mass):
    a_1 = 0.3; A = 1.785;
    return A * np.power(mass, 1.-a_1) / (1.-a_1)

def _cdf_midmass_Kroupa(mass):
    m_0 = 0.08; a_2 = 1.3; A = 1.785; B = 0.4352; k_1 = 0.08;
    return B + A * k_1 / (1.-a_2) * (np.power(mass, 1.-a_2) - np.power(m_0, 1.-a_2))
    
def _cdf_highmass_Kroupa(mass):
    m_1 = 0.5; a_3 = 2.3; A = 1.785; C = 0.86469; k_2 = 0.04;
    return C + A * k_2 / (1.-a_3) * (np.power(mass, 1.-a_3) - np.power(m_1, 1.-a_3))

# Takes a 2 element array representing the mass range for the provided IMF 
# distribution and produces n_samples number of stars in the given mass range
def generate_masses(mass_range,distrub='Chabrier',n_samples=10000):
    rang = np.arange(mass_range[0],mass_range[1],0.001)
    if distrub == 'Chabrier':
        m_0 = 1.;
        cdf = _cdf_lowmass_Chabrier(rang[rang <= m_0])
        cdf = np.append(cdf, _cdf_highmass_Chabrier(rang[rang > m_0]))
        inv_cdf = interpolate.interp1d(cdf, rang)
    elif distrub == 'Kroupa':
        m_0 = 0.08; m_1 = 0.5;
        cdf = _cdf_lowmass_Kroupa(rang[rang <= m_0])
        cdf = np.append(cdf, _cdf_midmass_Kroupa(rang[(rang > m_0) & (rang <= m_1)]))
        cdf = np.append(cdf, _cdf_highmass_Kroupa(rang[rang > m_1]))
        inv_cdf = interpolate.interp1d(cdf, rang)
    else:
        raise NameError("The " + distrub + " IMF is not provided in this method")
        return None 
    r = np.random.uniform(np.min(cdf),np.max(cdf),n_samples) 
    
    return inv_cdf(r)
    
def qdist(nsamples):
    #From Allen 2007
    gamma = 1.8
    randomq = powerlaw.rvs(gamma+1,size=nsamples)

    return randomq

