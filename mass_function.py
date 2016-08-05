# -*- coding: utf-8 -*-
"""
Created on Thu Aug 04 16:22:33 2016

@author: Caleb
"""

import pandas as pd
import numpy as np
from scipy import interpolate
#from astropy.table import Table
from scipy.special import erfc
import matplotlib.pyplot as plt
   
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
    return A * np.power(mass, 1-a_1) / (1 - a_1)

def _cdf_midmass_Kroupa(mass):
    m_0 = 0.08; a_2 = 1.3; A = 1.785; B = 0.4352; k_1 = 0.08;
    return B + A * k_1 / (1-a_2) * (np.power(mass, 1-a_2) - np.power(m_0, 1-a_2))
    
def _cdf_highmass_Kroupa(mass):
    m_1 = 0.5; a_3 = 2.3; A = 1.785; C = 0.86469; k_2 = 0.04;
    return C + A * k_2 / (1-a_3) * (np.power(mass, 1-a_3) - np.power(m_1, 1-a_3))

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
        print "The " + distrub + " IMF is not provided in this method"
        return None 
    r = np.random.uniform(np.min(cdf),np.max(cdf),n_samples) 
    return inv_cdf(r)