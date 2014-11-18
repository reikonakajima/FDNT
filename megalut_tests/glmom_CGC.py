"""
A demo that runs on the CGC branch ## (fields 005-9)

.. note:: this demo is valid for branches C**, R**.
"""

###################################################################################################
# Imports
#import megalut.great3
import megalut.meas
import megalut.sim
#import megalut.learn
import megalut.tools
import os.path

# Optional: set the logging level. If omitted, only warnings (and worse) will be shown.
import logging
logging.basicConfig(level=logging.INFO)
import numpy as np


###################################################################################################
# User-defined functions and classes needed for GREAT3
"""
def measure(img_fname,input_cat,stampsize,prefix):
	img = megalut.tools.image.loadimg(img_fname)
	meas_cat = megalut.meas.galsim_adamom.measure(img, input_cat, stampsize=stampsize,prefix=prefix)
	return meas_cat
"""

class CGC_simparams(megalut.sim.params.Params):
	def getrad(self):
		return np.random.uniform(0.7, 2.7)
	
	def getflux(self):
		if np.random.uniform() < 0.25:
			return np.random.uniform(70.0, 220.0)
		else:
			return np.random.uniform(10.0, 70.0)
	
	def getg(self):
		theta = np.random.uniform(0.0, 2.0* np.pi) 
		eps = np.random.uniform(0.0, 0.7) 
		return (eps * np.cos(2.0 * theta), eps * np.sin(2.0 * theta))
		
	def getsersicn(self, ix=0, iy=0, n=1):
		return 1.0 + (float(iy)/float(n)) * 2.0	
		# Lower sersic index = broader

"""
learnparams = megalut.learn.MLParams(
		name = "demo",
		features = ["gs_g1", "gs_g2", "gs_flux"],
		labels = ["tru_g1","tru_g2"],
		predlabels = ["pre_g1","pre_g2"],
		)

fannparams=megalut.learn.fannwrapper.FANNParams(
		hidden_nodes = [20, 20],
		max_iterations = 500,
	)
"""

basedir = "/vol/euclid2/euclid2_1/reiko/MegaLUT/temp/"

# Step 1: drawing the sims

simdir = os.path.join(basedir, "simdir")
simparams = megalut.sim.params.Params("CGC_simparams")
drawcatkwargs = {"n":30, "stampsize":64}
drawimgkwargs = {}

#megalut.sim.run.multi(simdir, simparams, drawcatkwargs, drawimgkwargs, ncat=2, nrea=3, ncpu=6)


# Step 2, measuring

import megalut.meas.glmomfunc

measdir = os.path.join(basedir, "meas_glmom")
measfct = megalut.meas.glmomfunc.measure
measfctkwargs = {}

#megalut.meas.run.onsims(simdir, simparams, measdir, measfct, measfctkwargs, ncpu=24, skipdone=False)


# Step 3, plotting comparison

import plot_diagnosis as pd

plot_pairs = [('tru_g1', 'glmom_g1'), ('tru_g2', 'glmom_g2'), ('tru_rad', 'glmom_sigma'),
	      ('tru_sersicn', 'glmom_rho4'), ('tru_flux', 'glmom_flux')]
pd.plot_diagnosis(simparams, measdir, plot_pairs)



# x, y, tru_g1, tru_g2, tru_rad, tru_sersicn, tru_sig,
# x, y, g1, g2, sigma, rho4, snratio



"""
###################################################################################################
# Start of the code

# Create an instance of the GREAT3 class
cgv=megalut.great3.great3.Run("control", "ground", "variable",
	datadir="/home/kuntzer/workspace/MegaLUT/great3_data_part",
	subfields=range(5,7))

# Now run the measurements on input images
cgv.meas("obs",measure,method_prefix="gs_")

# Make sim catalogs & images
cgv.sim(CGV_simparams(),n=10)

# Measure the observations with the same methods than the observation
cgv.meas("sim",measure,method_prefix="gs_")

# Train the ML
cgv.learn(learnparams=learnparams, mlparams=fannparams, method_prefix="gs_")

# Predict the output
cgv.predict(method_prefix="gs_")

# Write the output catalog
cgv.writeout("ML_FANN_demo_default")

# Prepare the presubmission file
# (This will fail as we work only on a subset of the data)
cgv.presubmit()
"""
