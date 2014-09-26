"""
A minimal demo about drawing simulated galaxies
"""

import logging
logging.basicConfig(level=logging.INFO)

import megalut
import megalut.sim
import sys
import fdnt_glmom
import numpy as np
import fdnt
import galsim

#
# debug flags
#
create_new_catalog_and_image = True


# First, we set the desired distributions of parameters, by overwriting the default distributions.

class MySimParams(megalut.sim.params.Params):
	def get_flux(self):
		return 300.0
	def get_g(self):
		return np.random.uniform(low=-0.05, high=0.05, size=2)
		
mysimparams = MySimParams()

if create_new_catalog_and_image:
	simcat = megalut.sim.stampgrid.drawcat(mysimparams, n=10, stampsize=48)
	megalut.utils.writepickle(simcat, "simcat.pkl")
else:
	simcat = megalut.utils.readpickle("simcat.pkl")

print simcat[:5]

# Now, we pass this catalog to drawimg, to generate the actual simulated images.

if create_new_catalog_and_image:
	megalut.sim.stampgrid.drawimg(simcat,
				      simgalimgfilepath="simgalimg.fits",
				      simtrugalimgfilepath="simtrugalimg.fits",
				      simpsfimgfilepath="simpsfimg.fits"
				      )

# We can directly proceed by measuring the images

gridimg = fdnt_glmom.loadimg("simpsfimg.fits")

meascat = fdnt_glmom.measure(gridimg, simcat, stampsize=48, prefix="mes")

# meascat is the output catalog, it contains the measured features:
print meascat[:5]

# We save it into a pickle
megalut.utils.writepickle(meascat, "meascat.pkl")

# Let's make a simple comparision plot:
import matplotlib.pyplot as plt
resi_x = meascat["mes_x"] - meascat["x"]
resi_y = meascat["mes_y"] - meascat["y"]
flag = meascat["mes_flag"]
plt.scatter(resi_x, resi_y, c=flag, lw=0, s=30)
plt.xlabel("mes_x residual")
plt.ylabel("mes_y residual")
plt.show()

