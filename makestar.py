def makestar(mass, age, model='Burrows97'):
    
    '''Calculates stellar properties such as Teff, radii, logg and logL using evolutionary models
    '''
    import pandas as pd
    import numpy as np
    from scipy import interpolate
    from astropy.table import Table
    import warnings

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
        burrows = pd.read_pickle("/Users/daniella/Python/Thesis/simulations/burrows97.pickle")
        allages = burrows["Age (Gyr)"]
        allmasses = burrows["M/Ms"]
        teff = burrows["Teff"]
        radius = burrows["R/Rs"]
        logg = burrows["logg(cgs)"]
        logL = burrows["logL/Ls"]
        
    if model == 'Baraffe03':
        #0.0005 - 0.1 Msun
        baraffe = pd.read_pickle("/Users/daniella/Python/Thesis/simulations/baraffe03.pickle")
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

#    starparams = [newteff,newrad,newlogg,newlogL]
    
    stardict = {'Teff (K)':newteff,'Radius (Rs)':newrad, 'log g':newlogg, 'log L':newlogL}    
    
    startable = Table(stardict)
    
    return startable
    
ex = makestar([0.06,0.09],[5,8])

print(ex)
