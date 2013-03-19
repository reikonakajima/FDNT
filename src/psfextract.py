#!/usr/bin/python

# a cheap script that converts PSFEx .psf (FITS) formats models to the ASCII format models the PSFExModel class can read

import pyfits
import sys
hdulist = pyfits.open(sys.argv[1])
tbldata = hdulist[1].data
data = tbldata.field(0)[0]

keywords=['LOADED','ACCEPTED','CHI2','POLNAXIS','POLGRP1','POLNAME1','POLZERO1','POLSCAL1','POLGRP2','POLNAME2','POLZERO2','POLSCAL2','POLNGRP','POLDEG1','PSF_FWHM','PSF_SAMP','PSFNAXIS','PSFAXIS1','PSFAXIS2','PSFAXIS3']

for i in keywords:
  try:
   c=hdulist[1].header[i]
  except:
   continue
  print i, c


print 'DATA'

for d in data:
 for dd in d:
  for ddd in dd:
   print ddd

hdulist.close()
