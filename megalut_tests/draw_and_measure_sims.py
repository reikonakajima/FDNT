"""
A minimal demo about drawing simulated galaxies
"""

import logging
logging.basicConfig(level=logging.INFO)

import megalut
import megalut.sim
import stampgrid
import sys
import fdnt_glmom
import numpy as np
import fdnt
import galsim

#
# debug flags
#
create_new_catalog_and_image = True

#
# run conditions
#
n = 10  # 1d size of the sample postage stamp array (there will be n x n galaxies)


# First, we set the desired distributions of parameters, by overwriting the default distributions.

class MySimParams(megalut.sim.params.Params):
	def get_sig(self):
		return 1.0

	def get_rad(self):
		#return np.random.uniform(0.5, 2.0)
		return 0.9

	def get_flux(self):
		#return np.random.uniform(10.0, 200.0)
		return 50.

	def get_g(self):
		#return np.random.uniform(low=-0.4, high=0.4, size=2)
		return (0., 0.)

	def get_sersicn(self, ix=0, iy=0, n=1):
		"""
		This is a bit special: we do not draw sersic indices randomly, as changing it
		from stamp to stamp significantly slows down galsim ! That's why we need to know
		the stamp index.

		:param ix: x index of the galaxy, going from 0 to n-1
		:param iy: y index, idem
		:param n: n x n is the number of stamps

		"""
		pseudorand = float(iy)/float(n)
		return 0.5 + pseudorand * 0.0


		
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

gridimg = fdnt_glmom.loadimg("simgalimg.fits")

meascat = fdnt_glmom.measure(gridimg, simcat, stampsize=48, prefix="mes")

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

