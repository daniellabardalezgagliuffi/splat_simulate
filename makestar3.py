from astropy.table import Table
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import numpy as np
import sys, glob
from astropy.io import ascii

def makestar(mass, age, model='Burrows97', test=False):

    Ages   = np.empty([0])
    Masses = np.empty([0])
    Temps  = np.empty([0])
    Lumins = np.empty([0])
    Radii  = np.empty([0])

    # Check for numpy arrays
    if isinstance(mass, np.ndarray) is False:
        mass = np.array([mass]).flatten()
    if isinstance(age, np.ndarray) is False:
        age  = np.array([age]).flatten()

    # Check that age and mass are the same size
    if len(age) != len(mass):
        raise Exception('Age and Mass must have the same number of values')

    inputs = np.array([(a,b) for a,b in zip(age, mass)])

    """
    for i in ['0.001', '0.005', '0.010', '0.050', '0.100',
              '0.120', '0.500', '1.000', '5.000', '10.000']:

        data = Table.read('Burrows/b97_%s'%i, comment='#', format='ascii')
        print float(i)
        #print data

        ages   = np.zeros(len(data)) + float(i)
        Ages   = np.append(Ages, ages)
        Masses = np.append(Masses, data['col1'].data)
        Temps  = np.append(Temps, data['col2'].data)
        Lumins = np.append(Lumins, 10**data['col3'].data)
        Radii  = np.append(Radii, data['col5'].data)
    """

    if model == 'Burrows97':

        for file in glob.glob('Burrows2/*'):

            if file.split('.')[1][-1] == 'n': 
                continue
            
            mass = float('0.'+file.split('.')[1])

            data = ascii.read(file, data_start=41, delimiter=' ')
            time = data['col2']
            temp = data['col3']
            lum  = data['col4']
            rad  = data['col5']*1e9 / 6.955e10

            # Just take the values with time > 0
            index    = np.where(time == 0)
            if len(index[0]) != 0:
                newindex = range(index[0], len(time))
            else: 
                newindex = range(len(time))

            Ages     = np.append(Ages, time[newindex])
            Masses   = np.append(Masses, np.zeros(len(newindex))+mass)
            Temps    = np.append(Temps, temp[newindex])
            Lumins   = np.append(Lumins, lum[newindex])
            Radii    = np.append(Radii, rad[newindex])

    #Luminosity  = griddata(np.array([Ages, Masses]).T, Lumins, (age, mass), method='linear', rescale=True)
    ##Temperature = griddata(np.array([Ages, Masses]).T, Temps, (age, mass), method='linear', rescale=True)
    #Radius      = griddata(np.array([Ages, Masses]).T, Radii, (age, mass), method='linear', rescale=True)

    if test:
        return Ages, Masses, Temps

    else:

        Luminosity  = griddata(np.array([Ages, Masses]).T, Lumins, inputs, method='linear', rescale=True)
        Temperature = griddata(np.array([Ages, Masses]).T, Temps, inputs, method='linear', rescale=True)
        Radius      = griddata(np.array([Ages, Masses]).T, Radii, inputs, method='linear', rescale=True)

        results = Table([Luminosity, Temperature, Radius], 
                     names = ['Luminosity', 'Temperature', 'Radius'])

        return results


def test1():
    ######## 3D plot

    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm
    import matplotlib.pyplot as plt
    import numpy as np

    Results = makestar(1,1,test=True)

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    #ax.plot_trisurf(Ages, Masses, Temps, cmap=cm.jet, linewidth=0.2, rasterized=True)
    ax.plot_trisurf(Results[0], Results[1], Results[2], cmap=cm.jet, linewidth=0.2, rasterized=True)
    ax.set_xlabel('Age')
    ax.set_ylabel('Mass')
    ax.set_zlabel('Luminosity')
    ax.set_zlabel('Temp')

    plt.show()

#test()
#sys.exit()

for mass in [0.2, 0.1, 0.05, 0.01, 0.005]:

    Ages, Masses = np.logspace(-3, 1, 10000), np.zeros(10000)+mass

    #print np.array([Ages, Masses]).T.shape
    #print Temps.shape

    Results = makestar(Masses, Ages)

    #grid_z1 = griddata(np.array([Ages, Masses]).T, Lumins, (grid_x, grid_y), method='linear', rescale=True)
    #print grid_z1
    #print Ages
    #print Results['Luminosity']

    plt.plot(Ages, Results['Luminosity'])

plt.xscale('log')
plt.yscale('log')
plt.xlabel('Age (Gyr)')
plt.ylabel('Luminosity (L_solar)')
plt.show()

