# -*- coding: utf-8 -*-
"""
Created on Thu Jul 21 17:03:06 2016

@author: daniella
"""

def invertTeffSpTrel():

    '''Never to be used again, but as record of how Teff2spt coefficients were found.''''

    import numpy as np
    from numpy.polynomial.polynomial import polyval, polyder
    import splat
    import scipy    

    #Allowed SpT range for Filippazzo et al. 2015 relation: 6-29 = M6-T9 
    spt = np.arange(47)*0.5+6
    Tcoeff = [4.747e+03, -7.005e+02, 1.155e+02, -1.191e+01, 6.318e-01, -1.606e-02, 1.546e-04]
    Teff = polyval(spt, Tcoeff)
    #polyval takes smallest degree first
    
    plt.plot(spt,Teff,'b')
    plt.xlabel('SpT')
    plt.ylabel('Teff')
    
    plt.plot(Teff,spt,'b')
    plt.ylabel('SpT')
    plt.xlabel('Teff')
    plt.savefig('filippazzo_teff_spt.jpg')
    
    #now to invert it, use Teff as x, spt as y and fit it
    
    poly1 = np.polyfit(Teff,spt,1)
    poly2 = np.polyfit(Teff,spt,2)
    poly3 = np.polyfit(Teff,spt,3)
    poly4 = np.polyfit(Teff,spt,4)
    poly5 = np.polyfit(Teff,spt,5)
    poly6 = np.polyfit(Teff,spt,6)
    poly7 = np.polyfit(Teff,spt,7)
    #polyfit returns highest degree first
    
    testspt1 = np.poly1d(poly1)
    testspt2 = np.poly1d(poly2)
    testspt3 = np.poly1d(poly3)
    testspt4 = np.poly1d(poly4)
    testspt5 = np.poly1d(poly5)
    testspt6 = np.poly1d(poly6)
    testspt7 = np.poly1d(poly7)
    ### poly1d takes highest degree coefficient first    
    
    plt.plot(Teff,spt,'--r')
    plt.plot(Teff,testspt1(Teff),'b')
    plt.plot(Teff,testspt2(Teff),'g')
    plt.plot(Teff,testspt3(Teff),'r')
    plt.plot(Teff,testspt4(Teff),'.b')
    plt.plot(Teff,testspt5(Teff),'.g')
    plt.plot(Teff,testspt6(Teff),'.r')
    plt.plot(Teff,testspt7(Teff),'--b')
    plt.legend(['original','deg1','deg2','deg3','deg4','deg5','deg6','deg7'])
    plt.ylabel('SpT')
    plt.xlabel('Teff')
    plt.savefig('polytry_spt_teff.jpg')
    
    chi4 = scipy.stats.chisquare(spt,testspt4(Teff))
    chi5 = scipy.stats.chisquare(spt,testspt5(Teff))
    chi6 = scipy.stats.chisquare(spt,testspt5(Teff))
    chi7 = scipy.stats.chisquare(spt,testspt5(Teff))
    
    #chi4 is the closest to 1    
    newcoeffs = poly4
    
    return newcoeffs
    

def spt2teff(spt):
    
    #Filippazzo et al. 2015 SpT to Teff relation
    #Valid for M6-T9
    import numpy as np
    from numpy.polynomial.polynomial import polyval, polyder
    import splat
    import scipy
    
    Tcoeff = [4.747e+03, -7.005e+02, 1.155e+02, -1.191e+01, 6.318e-01, -1.606e-02, 1.546e-04]
    Trms = 113
    Teff = polyval(spt, Tcoeff)   
    
    return teff
    

def teff2spt(teff):
    
    '''Calculates optical spectral type using Teff and Filippazzo et al. 2015 relation'''
    #based on filippazzo et al. 2015    
    newcoeffs = [ -5.38389120e-12, 3.73242761e-08, -8.74866546e-05, 6.82777082e-02, 1.20261522e+01]    
    sptfunc = np.poly1d(newcoeffs)
    sptn = sptfunc(teff)+10   #Burgasset SpTn notation
    
    spt = np.zeros(len(sptn))
    
    for i in range(len(sptn)):
        print(sptn[i])
        spt[i] = splat.typeToNum(sptn[i])
    
    return spt

plt.plot(spt,sptn-10,'.b')
plt.plot([5,30],[5,30],'r')
plt.xlabel('SpT (original)')
plt.ylabel('SpT (fit result)')
plt.savefig('sptvsspt.jpg')
