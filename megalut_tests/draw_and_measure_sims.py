"""
A minimal demo about drawing simulated galaxies
"""

import logging
logging.basicConfig(level=logging.INFO)

import megalut
import megalut.sim
import stampgrid
import sys
#import fdnt_glmom
import glmomfunc
import numpy as np
import fdnt
import galsim

#
# debug flags
#
create_new_catalog_and_image = False

#
# run conditions
#
n = 100  # 1d size of the sample postage stamp array (there will be n x n galaxies)


# First, we set the desired distributions of parameters, by overwriting the default distributions.

class MySimParams(megalut.sim.params.Params):

	## try a set of parameters similar to GREAT3 CGV parameters
        def get_rad(self):
                return np.random.uniform(0.7, 2.7)

        def get_flux(self):
                if np.random.uniform() < 0.25:
                        return np.random.uniform(70.0, 220.0)
                else:
                        return np.random.uniform(10.0, 70.0)

        def get_g(self):
                theta = np.random.uniform(0.0, 2.0* np.pi)
                eps = np.random.uniform(0.0, 0.7)
                return (eps * np.cos(2.0 * theta), eps * np.sin(2.0 * theta))

        def get_sersicn(self, ix=0, iy=0, n=1):
                return 1.0 + (float(iy)/float(n)) * 2.0
                # Lower sersic index = broader


		
mysimparams = MySimParams()

if create_new_catalog_and_image:
	simcat = stampgrid.drawcat(mysimparams, n=n, stampsize=48)
	megalut.tools.io.writepickle(simcat, "simcat.pkl")
else:
	simcat = megalut.tools.io.readpickle("simcat.pkl")

print simcat

# Now, we pass this catalog to drawimg, to generate the actual simulated images.

if create_new_catalog_and_image:
	stampgrid.drawimg(simcat,
			  simgalimgfilepath="simgalimg.fits",
			  simtrugalimgfilepath="simtrugalimg.fits",
			  simpsfimgfilepath="simpsfimg.fits"
			  )

# We can directly proceed by measuring the images

gridimg = megalut.tools.image.loadimg("simgalimg.fits")

meascat = glmomfunc.measure(gridimg, simcat, prefix="mes_")

# meascat is the output catalog, it contains the measured features:
print meascat

# We save it into a pickle
megalut.tools.io.writepickle(meascat, "meascat.pkl")

# Let's make a simple comparision plot:
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

resi_x = meascat["mes_x"] - meascat["x"]
resi_y = meascat["mes_y"] - meascat["y"]
flag = meascat["mes_flag"]
plt.scatter(resi_x, resi_y, c=flag, lw=0, s=30)
plt.xlabel("mes_x residual")
plt.ylabel("mes_y residual")
plt.show()

