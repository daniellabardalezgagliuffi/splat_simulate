import numpy as np

def ages(numstars, agerange=[0,10]):

	low, high = agerange
	ages = np.random.uniform(low, high, size=numstars)

	return ages
