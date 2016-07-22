import os
import numpy as np
from astropy.table import Table
import scipy.interpolate as interpolate

def densities(lffile = None, colors = None, interp = True, size = 10000):

    '''
    Reads in a file with a luminosity function and returns stellar density
    '''

    if lffile is None:
        lffile = os.path.dirname(os.path.realpath(__file__)) + \
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
