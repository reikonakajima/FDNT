#!/usr/bin/env python

import subprocess as S
import numpy as np
import pylab as P
import sys
import os
import random
import math
sys.path.append('/Users/reiko/2code/python/mylib')

"""
make_annz_input.py
  generates inputs for ANN (artificial neural network) code "ANNz" to generate output

"""

usage = """usage: %s <directory_name>
Generates inputs for the ANN code ANNz, for use with shear bias calibration.
<directory_name>: directory which contains the bias information for a given psf/gal/shear from 
                  the ring test "FDNTRing".  

                  The directory names are, e.g.,
                  moffatR15.exp.s01s00.p00p01/  (psf type+size, gal type, shear, psf shape),
                  where PSF is Moffat with ee50 of 1.5 pix, shear d1=0.01, and PSF shape e2=0.1.

                  The filenames in this directory have the format, e.g.,
                  moffatR15_e1_00_e2_01.sn80.exp_e1_09.dat (psf, sn, galaxy info),
                  with Moffat PSF of ee50=1.5 pix, S/N=80, and exponential galaxy with |e|=0.9.

output file #1:   The output file is "directory_name.dat", which compiles all data in the directory
                  into a single file.

output file #2,3,4: The output files are "directory_name.trn", "directory_name.vrf", and
                  "directory_name.tst", which takes the above file and puts it into
                  the ANNz input format of [data columns], [data error columns], [answers]
""" % (sys.argv[0]) 

# assumes that shear has been applied in the e1 direction only (magnitude truth_e1)
# and PSF is elliptical in e2 direction only (magnitude psf_e2)



def get_file_info_header_string(line, dir_name):
    """read file info from file header, as well as directory name
    """
    # get file info
    if line[0] != '#':
        raise IOError(f+' has incorrect data format')
    ls = line.split()
    psf_type = ls[1].split('_')[0].split('R')[0]
    psf_size = float(ls[1].split('_')[0].split('R')[1]) / 10.  # ee50 in pixel units
    psf_e1 = float(ls[1].split('_')[2]) / 10.
    psf_e2 = float(ls[1].split('_')[4]) / 10.
    sn = float(ls[2])
    gal_type = ls[3]
    try:
        gal_index = float(gal_type) / 10.
    except ValueError:
        if gal_type == 'exp':
            gal_index = 1.0
        elif gal_type == 'dev':
            gal_index = 4.0
        else:
            raise ValueError('unknown galaxy type '+gal_type)
    gal_e = float(ls[4])
    # get applied shear info from directory name
    dist1 = float(dir_name.split('.')[2].split('s')[1]) / 100.
    dist2 = float(dir_name.split('.')[2].split('s')[2]) / 100.

    return_str =  '%.2f %.2f %.2f  %.1f  %.2f %.2f  %.1f %.2f' \
        % (psf_size, psf_e1, psf_e2,  sn,  dist1, dist2,  gal_index, gal_e,)

    return return_str



if len(sys.argv) != 2 or not os.path.isdir(sys.argv[1]):
    print usage
    sys.exit(1)

dir_name = sys.argv[1]
out_filename = dir_name + '.dat'
outf = open(out_filename, 'w')

shcommand = 'ls %s/*[0-9].dat' % (dir_name)   # this is the bias data
file_list = S.Popen(shcommand, shell=True, stdout=S.PIPE).communicate()[0].split()


for f in file_list:
    is_first_line = True
    for line in open(f).readlines():
        if is_first_line:
            header_str = get_file_info_header_string(line, dir_name)
            is_first_line = False
            ls = line.split()
        else:
            if line[0] == '#':  # indicates missing data
                continue
            print >> outf, header_str, line.strip()

outf.close()

truth_e1 = float(dir_name.split('.')[2].split('s')[1]) / 100.
truth_e2 = float(dir_name.split('.')[2].split('s')[2]) / 100.
psf_e1 = float(ls[1].split('_')[2]) / 10.
psf_e2 = float(ls[1].split('_')[4]) / 10.

print truth_e1,truth_e2
print psf_e1,psf_e2

if(psf_e1 != 0 or truth_e2 != 0 or truth_e1==0):
 print "I cannot work with theses parameters"
 sys.exit(1)


train_filename = out_filename.replace('.dat', '.train')
print "train filename",train_filename
trainf = open(train_filename, 'w')
#verification_filename = out_filename.replace('.dat', '.verify')
#verificationf = open(verification_filename, 'w')
#test_filename = out_filename.replace('.dat', '.test')
#testf = open(test_filename, 'w')

outf = {(1,2,3,4,5,6): trainf } #, 
 #       (): verificationf, 
 #       (): testf,
 #       }

err_str = '0 0 0'   # errors on S/N, gal |e|, and gal ee50, not needed for training
for line in open(out_filename).readlines():
    
    (p1, p2, p3, sn, d1, d2, gi, gal_e, gal_size, gal_ee50, \
         mean1, mean2, errm1, errm2, r1, r2, n1, n2) = line.split()
    
    randindex = random.randint(1,6)
    for key in outf.keys():
        if randindex in key:
            break
    if(float(gal_ee50)>1.5 and float(gal_e)<0.85):
      print >> outf[key],  1./float(sn)/float(sn), gal_e, math.exp(-float(gal_ee50)),  err_str,  float(mean1)/truth_e1, float(errm1)/truth_e1


trainf.close()
#verificationf.close()
#testf.close()
