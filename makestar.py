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

    newteff = interpteff(age,mass)
    newrad = interprad(age,mass)
    newlogg = interplogg(age,mass)
    newlogL = interplogL(age,mass)

    starparams = [newteff,newrad,newlogg,newlogL]
    
    stardict = {'Teff (K)':newteff,'Radius (Rs)':newrad, 'log g':newlogg, 'log L':newlogL}    
    
    startable = Table(stardict)
    
    return startable
    
ex = makestar(5,0.06)

print(ex)
