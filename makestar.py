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
