import pandas as pd
import numpy as np
from scipy import interpolate
from astropy.table import Table

def makestar(age, mass, model='Burrows97'):
    
    if model == 'Burrows97':

        burrows = pd.read_pickle("burrows97.pickle")
        allages = np.array(burrows["Age (Gyr)"]) 
        allmasses = np.array(burrows["M/Ms"]) 
        teff = np.array(burrows["Teff"])
        radius = np.array(burrows["R/Rs"])
        logg = np.array(burrows["logg(cgs)"])
        logL = np.array(burrows["logL/Ls"])
        
    if model == 'Baraffe03':
        
        baraffe = pd.read_pickle("baraffe03.pickle")
        allages = np.array(baraffe["Age (Gyr)"]) 
        allmasses = np.array(baraffe["M/Ms"]) 
        teff = np.array(baraffe["Teff"])
        radius = np.array(baraffe["R/Rs"])
        logg = np.array(baraffe["logg(cgs)"])
        logL = np.array(baraffe["logL/Ls"])

    interpteff = interpolate.interp2d(allages,allmasses,teff,kind='linear')
    interprad = interpolate.interp2d(allages,allmasses,radius,kind='linear')
    interplogg = interpolate.interp2d(allages,allmasses,logg,kind='linear')
    interplogL = interpolate.interp2d(allages,allmasses,logL,kind='linear')

    if len(age) > 1:
        newteff = np.zeros(len(age))
        newrad = np.zeros(len(age))
        newlogg = np.zeros(len(age))
        newlogL = np.zeros(len(age))
        for i in range(len(age)):
            newteff[i] = interpteff(age[i],mass[i])
            newrad[i] = interprad(age[i],mass[i])
            newlogg[i] = interplogg(age[i],mass[i])
            newlogL[i] = interplogL(age[i],mass[i])
    else:
        newteff = interpteff(age,mass)
        newrad = interprad(age,mass)
        newlogg = interplogg(age,mass)
        newlogL = interplogL(age,mass)
        
    stardict = {'Mass (Ms)':mass, 'Age (Gyr)':age, 'Teff (K)':newteff,'Radius (Rs)':newrad, 'log g':newlogg, 'log L':newlogL}
    
    startable = Table(stardict)
    
    return startable
    
    
# Chabrier CDF taken from Jumper & Fisher (2013) w/ params from Chabrier (2005)
def _cdf_lowmass_Chabrier(mass):
    A_s = 0.724;  sig = 0.55; m_c = 0.2;
    return A_s * sig * np.sqrt(np.pi/2.) * erfc((np.log10(m_c) - np.log10(mass)) / (sig * np.sqrt(2.0)))
def _cdf_highmass_Chabrier(mass):
    B_s = 0.323; x = 1.35; m_o = 1.; C = 0.896;
    return C + (B_s / np.log(10.)) * (np.power(m_o,-x) / x) * (1. - np.power((mass / m_o),-x))
    
def _cdf_lowmass_Kroupa(mass):
    a_1 = 0.3; A = 0.194;
    return A * np.power(mass, 1-a_1) / (1 - a_1)

def _cdf_midmass_Kroupa(mass):
    m_0 = 0.08; a_2 = 1.3; A = 0.194; B = 0.0474;
    return B + A / (1-a_2) * (np.power(mass, 1-a_2) - np.power(m_0, 1-a_2))
    
def _cdf_highmass_Kroupa(mass):
    m_1 = 0.5; a_3 = 2.3; A = 0.194; C = 0.632; 
    return C + A / (1-a_3) * (np.power(mass, 1-a_3) - np.power(m_1, 1-a_3))

def generate_masses(mass_range,distrub='Chabrier',n_samples=10000):
    rang = np.arange(mass_range[0],mass_range[1],0.001)
    if distrub == 'Chabrier':
        m_0 = 1.;
        cdf = _cdf_lowmass_Chabrier(rang[rang <= m_0])
        cdf = np.append(cdf, _cdf_highmass_Chabrier(rang[rang > m_0]))
        inv_cdf = interpolate.interp1d(cdf, rang)
    # This one may be slightly off
    elif distrub == 'Kroupa':
        m_0 = 0.08; m_1 = 0.5;
        cdf = _cdf_lowmass_Kroupa(rang[rang <= m_0])
        cdf = np.append(cdf, _cdf_midmass_Kroupa(rang[(rang > m_0) & (rang <= m_1)]))
        cdf = np.append(cdf, _cdf_highmass_Kroupa(rang[rang > m_1]))
        inv_cdf = interpolate.interp1d(cdf, rang)
    else:
        print "The " + distrub + " IMF is not provided in this method"
        return None 
    r = np.random.uniform(np.min(cdf),np.max(cdf),n_samples) 
    return inv_cdf(r)

