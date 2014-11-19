"""
A demo that tests FDNT with megalut sims and measurements

"""

###################################################################################################
# Imports
import megalut.meas
import megalut.sim
import megalut.tools
import os.path

# Optional: set the logging level. If omitted, only warnings (and worse) will be shown.
import logging
logging.basicConfig(level=logging.INFO)
import numpy as np


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

basedir = "/vol/euclid2/euclid2_1/reiko/MegaLUT/temp/"

# Step 1: drawing the sims

simdir = os.path.join(basedir, "simdir")
simparams = CGC_simparams()
drawcatkwargs = {"n":30, "stampsize":64}
drawimgkwargs = {}

#megalut.sim.run.multi(simdir, simparams, drawcatkwargs, drawimgkwargs, ncat=2, nrea=3, ncpu=6, savepsfimg=True, savetrugalimg=True)


# Step 2, measuring

import megalut.meas.fdntfunc

measdir = os.path.join(basedir, "meas_fdnt")
measfct = megalut.meas.fdntfunc.measure
measfctkwargs = {"sewpy_workdir": os.path.join(measdir,'sewpy'), "stampsize": 64}

megalut.meas.run.onsims(simdir, simparams, measdir, measfct, measfctkwargs, ncpu=1, skipdone=False)


# Step 3, plotting comparison

import plot_diagnosis as pd

plot_pairs = [('tru_g1', 'fdnt_g1'), ('tru_g2', 'fdnt_g2'), ('tru_rad', 'fdnt_sigma'),
	      ('tru_sersicn', 'fdnt_b22'), ('tru_flux', 'fdnt_flux')]
pd.plot_diagnosis(simparams, measdir, plot_pairs)



# x, y, tru_g1, tru_g2, tru_rad, tru_sersicn, tru_sig,
# x, y, g1, g2, sigma, rho4, snratio
