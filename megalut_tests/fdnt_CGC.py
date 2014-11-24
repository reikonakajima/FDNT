"""
A demo that runs on the CGC branch ## (fields 005-9)

.. note:: this demo is valid for branches C**, R**.
"""

###################################################################################################
# Imports
import megalut.great3
import megalut.meas
import megalut.sim
import megalut.learn
import megalut.tools
import os.path
import megalut.meas.fdntfunc

# Optional: set the logging level. If omitted, only warnings (and worse) will be shown.
import logging
logging.basicConfig(level=logging.INFO)
import numpy as np

from megalut.meas import fdnt_measfct

###################################################################################################
# User-defined functions and classes needed for GREAT3

class CGC_simparams(megalut.sim.params.Params):
	def get_rad(self):
		return np.random.uniform(0.7, 2.7)
	
	def get_flux(self):
		if np.random.uniform() < 0.25:
			return np.random.uniform(70.0, 220.0)
		else:
			return np.random.uniform(10.0, 70.0)

	def get_g(self):
		theta = np.random.uniform(0.0, 2.0*np.pi)
		eps = np.random.uniform(0.0, 0.7) 
		return (eps*np.cos(2.0*theta), eps*np.sin(2.0*theta))
		
	def get_sersicn(self, ix=0, iy=0, n=1):
		return 1.0 + (float(iy)/float(n)) * 2.0	
		# Lower sersic index = broader

learnparams = megalut.learn.MLParams(
		name = "demo",
		features = ["fdnt_g1", "fdnt_g2", "fdnt_flux"],
		labels = ["tru_g1","tru_g2"],
		predlabels = ["pre_g1","pre_g2"],
		)

fannparams=megalut.learn.fannwrapper.FANNParams(
		hidden_nodes = [20, 20],
		max_iterations = 500,
	)

basedir = "/vol/euclid2/euclid2_1/reiko/MegaLUT/temp/"

simparam_name = 'CGC_FDNT_debug'
cgc_simparam = CGC_simparams(simparam_name)

"""
import plot_diagnosis as pd
plot_pairs = [('tru_g1', 'fdnt_g1'), ('tru_g2', 'fdnt_g2'), ('tru_rad', 'fdnt_sigma'),
	      ('tru_sersicn', 'fdnt_b22'), ('tru_flux', 'fdnt_flux')]
pd.plot_diagnosis(simparams, measdir, plot_pairs)
"""

###################################################################################################
# Start of the code

# Create an instance of the GREAT3 class
cgc = megalut.great3.great3.Run("control", "ground", "constant",
	datadir="/vol/euclid4/euclid4_1/reiko/GREAT3_DATA/",
	subfields=range(5,6))

# Now run the measurements on input images
cgc.meas_obs(fdnt_measfct.measfct, skipdone=False, ncpu=1)

exit()

"""
# Make sim catalogs & images
cgc.sim(CGC_simparams(), n=10)

# Measure the observations with the same methods than the observation
cgc.meas("sim", measure, method_prefix="fdnt_")

# Train the ML
cgc.learn(learnparams=learnparams, mlparams=fannparams, method_prefix="gs_")

# Predict the output
cgc.predict(method_prefix="gs_")

# Write the output catalog
cgc.writeout("ML_FANN_demo_default")

# Prepare the presubmission file
# (This will fail as we work only on a subset of the data)
cgc.presubmit()
"""
