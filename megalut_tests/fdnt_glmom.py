"""
Shape measurement with FDNT's GLMoments.
"""

import numpy as np
import sys, os
import copy
from datetime import datetime

import logging
logger = logging.getLogger(__name__)

import astropy.table
import galsim
import fdnt



def loadimg(img_file_path):
	"""
	Load an image from a FITS file into GalSim::Image, enforcing that the origin is (0, 0).
	
	:param imgfilepath: path to FITS image
	:returns: GalSim Image
	"""
	
	logger.info("Loading FITS image %s..." % (os.path.basename(img_file_path)))
	bigimg = galsim.fits.read(img_file_path)
	# note: no need to setOrigin() here for FDNT (postage stamp extraction not needed)
	logger.info("Done with loading %s, shape is %s" % (os.path.basename(img_file_path),
							   bigimg.array.shape))
	
	logger.warning("The origin and stampsize conventions are new and should be tested !")
	
	return bigimg



def measure(bigimg, catalog, xname="x", yname="y", stampsize=100, prefix="mes_glmom"):
	"""
	Use the pixel positions provided via the input table to extract postage stamps from
	the image and measure their shape parameters.
	Return a copy of your input catalog with the new columns appended.
	One of these columns is the flag:
	
	* 0: OK
	* 1: stamp is not fully within image  (not implemented for GLMoments)
	* 2: galsim centroid failed (distance > 10 pixels...)
	* 3: galsim failed
	
	:param bigimg: galsim (postage stamp mosaic) image
	:param catalog: astropy table of objects to be measured
	:param xname: column name containing the x coordinates in pixels
	:param yname: column name containing the y coordinates in pixels
	:param stampsize: width = height of stamps, must be even
	:type stampsize: int
	:param prefix: a string to prefix the field names that I'll write
		
	:returns: astropy table
	
	"""
	
	assert int(stampsize)%2 == 0 # checking that it's even
	
	starttime = datetime.now()
	logger.info("Starting shape measurements on %ix%i stamps" % (stampsize, stampsize))
	
	# Prepare an output table with all the required columns
	output = copy.deepcopy(catalog)
	output.add_columns([
		astropy.table.Column(name=prefix+"_flag", data=np.zeros(len(output), dtype=int)),
		astropy.table.Column(name=prefix+"_gl_flag", data=np.zeros(len(output), dtype=int)),
		astropy.table.Column(name=prefix+"_flux", dtype=float, length=len(output)),
		astropy.table.Column(name=prefix+"_x", dtype=float, length=len(output)),
		astropy.table.Column(name=prefix+"_y", dtype=float, length=len(output)),
		astropy.table.Column(name=prefix+"_g1", dtype=float, length=len(output)),
		astropy.table.Column(name=prefix+"_g2", dtype=float, length=len(output)),
		astropy.table.Column(name=prefix+"_sigma", dtype=float, length=len(output)),
		astropy.table.Column(name=prefix+"_rho4", dtype=float, length=len(output)),
		astropy.table.Column(name=prefix+"_snratio", dtype=float, length=len(output)),
	])
	
	# Can have boolean columns:
	#astropy.table.Column(name=prefix+"_ok", data=np.zeros(len(output), dtype=bool))
	# Could use MaskedColumn here !
	#astropy.table.MaskedColumn(name="mes_adamom_flux", dtype=float, length=len(output),
	#			   fill_value=-1)
	
	# Save something useful to the meta dict
	#output.meta[prefix + "_imgfilepath"] = bigimg.origimgfilepath
	output.meta[prefix + "_xname"] = xname
	output.meta[prefix + "_yname"] = yname
	
	n = len(output)
	
	print "column names of output", output.colnames

	# Loop over each object
	for gal in output:
		
		# Some simplistic progress indication:
		if gal.index%5000 == 0:  # is "index" an astropy table entry?
			logger.info("%6.2f%% done (%i/%i) " % (100.0*float(gal.index)/float(n),
							       gal.index, n))
		# get centroid from catalog
		(x, y) = (gal[xname], gal[yname])
		g1g2 = (gal['tru_g1'], gal['tru_g2'])

		# We measure the moments... GLMoment may fail from time to time, hence the try:
		try:
			# TODO:  GLMoments tend to fail for sigma < 2.71828 (e).  FIX!!
			res = fdnt.GLMoments(bigimg, x, y, gal['tru_rad'], guess_g1g2=g1g2)
			
		except RuntimeError, m:
			print m
			# This is awesome, but clutters the output 
			#logger.exception("GLMoments failed on: %s" % (str(gal)))
			# So insted of logging this as an exception, we use debug, but include
			# the traceback :
			logger.debug("GLMoments failed with %s:\n %s" % (m, str(gal)), exc_info=True)
			#print "GLMoments failed on:\n %s" % (str(gal))
			gal[prefix + "_flag"] = 3
			continue

		try:
			s = galsim.Shear(e1=res.observed_e1, e2=res.observed_e2)
			g1 = s.getG1()
			g2 = s.getG2()
		except ValueError:
			g1 = res.observed_e1   # they should be "nonsense value" == -10.
			g2 = res.observed_e2

		gal[prefix+"_flux"] = res.observed_b00
		gal[prefix+"_x"] = res.observed_centroid.x
		gal[prefix+"_y"] = res.observed_centroid.y
		gal[prefix+"_g1"] = g1
		gal[prefix+"_g2"] = g2
		gal[prefix+"_sigma"] = res.observed_sigma
		gal[prefix + "_gl_flag"] = res.observed_flags
		# note: b_22 = rho4-4*rho2+2 = rho4-4*b_11+2*b_00;  b22/b00 is a substitute
		gal[prefix+"_rho4"] = res.observed_b22/res.observed_b00
		gal[prefix + "_snratio"] = res.observed_significance

		# If we made it so far, we check that the centroid is roughly ok:
		if np.hypot(x - gal[prefix+"_x"], y - gal[prefix+"_y"]) > 10.0:
			gal[prefix + "_flag"] = 2
		
	endtime = datetime.now()	
	logger.info("All done")

	nfailed = np.sum(output[prefix+"_flag"] > 0)
	
	logger.info("I failed on %i out of %i sources (%.1f percent)" % \
			    (nfailed, n, 100.0*float(nfailed)/float(n)))
	logger.info("This measurement took %.3f ms per galaxy" % \
			    (1e3*(endtime - starttime).total_seconds() / float(n)))
	
	return output



		
def getstamp(x, y, bigimg, stampsize):
	"""
	Prepare a bounded galsim image stamp "centered" at position (x, y) of the input galsim image.
	Can use the array attribute of the stamp if you want to get the actual pixels.
	
	:returns: a tuple(GalSim.Image, flag)
	"""

	assert int(stampsize)%2 == 0 # checking that it's even

	xmin = int(np.round(x-0.5)) - int(stampsize)/2
	xmax = int(np.round(x-0.5)) + int(stampsize)/2 - 1
	ymin = int(np.round(y-0.5)) - int(stampsize)/2
	ymax = int(np.round(y-0.5)) + int(stampsize)/2 - 1
			
	assert ymax - ymin == stampsize - 1
	assert xmax - xmin == stampsize - 1
	
	# We check that these bounds are fully within the image
	if xmin < bigimg.getXMin() or xmax > bigimg.getXMax() \
		    or ymin < bigimg.getYMin() or ymax > bigimg.getYMax():
		return (None, 1)
		
	# We prepare the stamp
	bounds = galsim.BoundsI(xmin, xmax, ymin, ymax)
	stamp = bigimg[bounds] # galaxy postage stamp
	assert stamp.array.shape == (stampsize, stampsize)
	
	return (stamp, 0)
	
	
	
	
	
#def npstampgrid(bigimg, catalog, xname="x", yname="y", stampsize=100):
#	"""
#	I build a numpy array with stamps, intended for visualization
#	"""
#
#	#n = len(catalog)
#	#nrows = int(np.ceil(n/10))
#	#big = np.zeros((10*stampsize, nrows*stampsize))
#	#print n, nrows
#	
#	stamplist = []
#	for gal in catalog:
#		(x, y) = (gal[xname], gal[yname])
#		(gps, flag) = getstamp(x, y, bigimg, stampsize)
#		if flag != 0:
#			stamplist.append(np.zeros(stampsize, stampsize))
#		else:
#			stamplist.append(gps.array)
#	
#	big = np.vstack(stamplist)
#	return big



def pngstampgrid(pngfilepath, bigimg, catalog, xname="x", yname="y", stampsize=100, ncols=5,
		 upsample=4, z1="auto", z2="auto"):
	"""
	I write a grid of stamps corresponding to your catalog in a png image, so that
	you can visualize those galaxies...
	For this I need f2n, but anyway I'm just a little helper.
	"""
	
	try:
		import f2n
	except:
		logger.exception("Could not import f2n, will skip this...")
		return
		
	n = len(catalog)
	nrows = int(np.ceil(float(n)/float(ncols)))
		
	stamprows = []
	for nrow in range(nrows):
		stamprow = []
		for ncol in range(ncols):
			
			index = ncol + ncols*nrow
			if index < n: # Then we have a galaxy to show
				gal = catalog[index]
				(x, y) = (gal[xname], gal[yname])
				(gps, flag) = getstamp(x, y, bigimg, stampsize)
				npstamp = gps.array
				
				f2nstamp = f2n.f2nimage(numpyarray=npstamp, verbose=False)

				f2nstamp.setzscale(z1, z2)
				f2nstamp.makepilimage("log", negative = False)
				f2nstamp.upsample(upsample)
			
				txt = [
					"%i (%i, %i)" % (gal.index, x, y),
				]
				f2nstamp.writeinfo(txt, colour=(255, 255, 100))
				
			else: # No more galaxies, we just fill the splot with a grey empty stamp.
				npstamp = np.zeros((stampsize, stampsize))
				f2nstamp = f2n.f2nimage(numpyarray=npstamp, verbose=False)
				f2nstamp.setzscale(-1.0, 1.0)
				f2nstamp.makepilimage("lin", negative = False)
				f2nstamp.upsample(4)
			
			
			stamprow.append(f2nstamp)
		stamprows.append(stamprow)
	f2n.compose(stamprows, pngfilepath)
	logger.info("Wrote %s" % (pngfilepath))
