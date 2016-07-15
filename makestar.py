def makestar(age, mass, model):
    
    import pandas as pd
    import numpy as np
    from scipy import interpolate
    
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

    interpteff = scipy.interpolate.interp2d(allages,allmasses,teff,kind='linear')
    interprad = scipy.interpolate.interp2d(allages,allmasses,radius,kind='linear')
    interplogg = scipy.interpolate.interp2d(allages,allmasses,logg,kind='linear')
    interplogL = scipy.interpolate.interp2d(allages,allmasses,logL,kind='linear')

    newteff = interpteff(age,mass)
    newrad = interprad(age,mass)
    newlogg = interplogg(age,mass)
    newlogL = interplogL(age,mass)

    starparams = [newteff,newrad,newlogg,newlogL]
                
    return starparams
