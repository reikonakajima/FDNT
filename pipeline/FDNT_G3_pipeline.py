#!/usr/bin/env python

import subprocess as S
import numpy as np
import pylab as P
import sys
#sys.path.append('/Users/reiko/2code/python/mylib')

"""
FDNT_G3_pipeline.py
  GREAT3 pipeline for FDNT analysis

  Take first argument as resolution value
"""

if len(sys.argv) == 1:
    print 'first argument specifies a key to the file to be read'
    sys.exit(1)

#gal_num = sys.argv[1]
#gal_dilation = sys.argv[2]


# FITS image lists
#gal_fits_image_list = ['output/realgals_moffat_d%s_gt%s.fits' % (gal_dilation, gal_num),]
#gal_fits_image_list = ['output/g08_gal_R20.fits',]
gal_fits_image_list = ['%s' % (sys.argv[1],),]
#gal_fits_image_list = ['multi.fits',]
#psf_fits_image_list = ['g08_psf.fits',]

# SExtractor parameters
#se_config_file = '/home/reiko/2code/se/g3_high_sn.sex'
#se_output_param_file = '/home/reiko/2code/se/glfit.param'
se_config_file = '/users/marggraf/src/Euclid/FDNT/fdnt1/1_pipeline/pipeline/g3_high_sn.sex'
se_output_param_file = '/users/marggraf/src/Euclid/FDNT/fdnt1/1_pipeline/pipeline/glfit.param'
se_conv_file = None   # should depend on the given PSF



for i, gal_image in enumerate(gal_fits_image_list):

    #
    # set parameters
    #
    h = gal_image.split('/')[-1].split('.fits')[0]
    se_output_catalog = '%s.se.cat' % (h,)
    ds9_reg_se_fname = '%s.se.reg' % (h,)
    ds9_reg_fname = '%s.reg' % (h,)

    #psf_image = psf_fits_image_list[i]
    se_psf_output_catalog = '%s.psf.se.cat' % (h,)
    ds9_psf_reg_fname = '%s.psf.se.reg' % (h,)
    #psf_bvec_filename = '%s.psf.bvec' % (h,)
    psf_filename = 'moffat'  # this keyword triggers use of the known PSF of GREAT08

    print 'working on', h, '...'

    #
    # *** GET ORIGINAL SExtractor GALAXY CATALOG ***
    #
    se_options = ' -CATALOG_NAME %s' % (se_output_catalog,)
    se_options += ' -PARAMETERS_NAME %s' % (se_output_param_file,)
    se_options += ' -FILTER N'
    se_shcommand = '/scisoft/bin/sex %s -c %s %s' % (gal_image, se_config_file, se_options)
    S.call(se_shcommand, shell=True)
    #
    # get SE image regions
    # 
    shcommand = '/users/marggraf/src/Euclid/MakeDS9Ellipses/MakeDS9Ellipses_se < %s > %s' % (se_output_catalog, ds9_reg_se_fname,)
    S.call(shcommand, shell=True)

    #
    # *** GET PSF BVECTOR ***
    #
    # use known PSF info
    #
    """
    # Run SExtractor
    #
    se_options = ' -CATALOG_NAME %s' % (se_psf_output_catalog,)
    se_options += ' -PARAMETERS_NAME %s' % (se_output_param_file,)
    se_options += ' -FILTER N'
    se_shcommand = '/scisoft/bin/sex %s -c %s %s' % (psf_image, se_config_file, se_options)
    S.call(se_shcommand, shell=True)
    #
    # get SE image regions
    # 
    shcommand = '/home/reiko/bin/MakeDS9Ellipses_se < %s > %s' % (se_psf_output_catalog, ds9_psf_reg_fname,)
    S.call(shcommand, shell=True)
    #
    # fit 2d GL polynomials
    #
    shcommand = ''
    """

    #
    # *** RUN FDNT ***
    #
    # Generate FDNTGreat parameter input ASCII file
    #
    fdnt_param_fname = '%s.par' % (h,)
    stampSize = 80
    sky = 0.
    readNoise = 1000.  # flat noise!  setting this to zero will cause the code to puke.
    gain = 1.
    idCol       = 1                   # Object ID
    xCol        = 2                   # x centroid
    yCol        = 3                   # y centroid
    magCol      = 4                   # Magnitude
    rCol        = 8                   # Half-light radius
    aCol        = 9                   # Major axis
    bCol        = 10                  # Minor axis
    paCol       = 11                  # Position angle
    parfile_content =  'fitsName = %s[0]\n' % (gal_image,)
    parfile_content += 'catName = %s\n' % (se_output_catalog,)
    parfile_content += 'psfName = %s\n' % (psf_filename,)
    parfile_content += 'stampSize = %d\n' % (stampSize,)
    parfile_content += 'sky = %f\n' % (sky,)
    parfile_content += 'rn = %f\n' % (readNoise,)
    parfile_content += 'gain = %f\n' % (gain,)
    #parfile_content += ' = %s\n' % ()
    fpf = open(fdnt_param_fname, 'w')
    print >> fpf, parfile_content
    fpf.close()
    #
    # Process FDNT output
    #
    fdnt_catalog = '%s.out' % (h,)
#    fdnt_shcommand = '/home/reiko/bin/FDNTGreat_KillNyquist < %s > %s' % (fdnt_param_fname, fdnt_catalog)
    fdnt_shcommand = '/users/marggraf/src/Euclid/FDNT/fdnt1/1_pipeline/src/FDNTGreat < %s > %s' % (fdnt_param_fname, fdnt_catalog)
    S.call(fdnt_shcommand, shell=True)

    #
    # get SE image regions
    # 
    shcommand = 'MakeDS9Ellipses_FDNTGreat < %s > %s' % (fdnt_catalog, ds9_reg_fname,)
    S.call(shcommand, shell=True)
