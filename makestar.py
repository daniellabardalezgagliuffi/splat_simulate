import pandas as pd
import numpy as np
from scipy import interpolate
from astropy.table import Table
import warnings

def makestar(age, mass, model='Burrows97'):

    if np.min(mass) < 0.0005:
        raise NameError('Mass below minimum mass of 0.0005Msun')

    if np.max(mass) > 0.1 and model=='Baraffe03':
       warnings.warn('Mass above maximum mass of 0.1Msun for Baraffe 2003. Using Burrows 1997 instead.')  
 #       print('Mass above maximum mass of 0.1Msun for Baraffe 2003. Using Burrows 1997 instead.')
       model = 'Burrows97'

    if np.min(mass) > 0.2:
        raise NameError('Mass above maximum mass of 0.2Msun for Burrows 1997')    

    if model == 'Burrows97':
        #0.0005 - 0.2 Msun
        burrows = pd.read_pickle("burrows97.pickle")
        allages = np.array(burrows["Age (Gyr)"]) 
        allmasses = np.array(burrows["M/Ms"]) 
        teff = np.array(burrows["Teff"])
        radius = np.array(burrows["R/Rs"])
        logg = np.array(burrows["logg(cgs)"])
        logL = np.array(burrows["logL/Ls"])
        
    if model == 'Baraffe03':
        #0.0005 - 0.1 Msun
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
