def ages(numstars, agerange=[0,10]):

	import numpy as np

	low, high = agerange
	ages = np.random.uniform(low, high, size=numstars)

	return ages
