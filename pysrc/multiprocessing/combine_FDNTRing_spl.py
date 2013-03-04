#!/usr/bin/env python

import subprocess as S
import numpy as np
import pylab as P
import sys
sys.path.append('/Users/reiko/2code/python/mylib')
import shearestimator as se

"""
combine_FDNTRing_spl.py
  combine the FDNTRing_spl outputs

"""

# collect the set of outputs
header = 'test'
shcommand = 'ls %s[0-9]*.out' % (header)
dataflist = S.Popen(shcommand, shell=True, stdout=S.PIPE).communicate()[0].strip().split()

# check the set in datalist is correct
dataf = open(dataflist[0], 'r')
nTheta_ref = int(filter(lambda s: '#nTheta' in s, dataf.readlines())[0].split()[2])
print 'nTheta is', nTheta_ref
dataf.close()
theta_vals = np.zeros(nTheta_ref)
nMeasuredTheta = np.zeros(nTheta_ref)
nFailTheta = np.zeros(nTheta_ref)

# add up results
seMean = se.ShearEstimator()
seRaw = se.ShearEstimator()
seNoFix = se.ShearEstimator()
for datafname in dataflist:

    dataf = open(datafname, 'r')
    file_top = dataf.tell()
    nTheta = int(filter(lambda s: '#nTheta' in s, dataf.readlines())[0].split()[2])
    if nTheta != nTheta_ref:
        print 'new nTheta is', nTheta
        print 'data file %s incompatible with nTheta = %d' % (datafname, nTheta_ref)
        sys.exit(1)
    dataf.seek(file_top)
    data = map(lambda s: s.strip(), filter(lambda s: s[0] != '#', dataf.readlines()))

    counts_data = data[0].split()
    iTheta = int(counts_data[0])
    theta_vals[iTheta] = float(counts_data[1])
    nMeasuredTheta[iTheta] += int(counts_data[2])
    nFailTheta[iTheta] += int(counts_data[3])

    seMean_data = data[1]
    sum1, sum2, sumsq1, sumsq2, N, sumFixEst, theta = [float(f) for f in seMean_data.split()]
    mean1 = sum1 / N
    mean2 = sum2 / N
    temp1 = sumsq1/N - mean1*mean1
    if temp1 <= 0.:
        sig1 = 0.
    else:
        sig1 = np.sqrt(temp1)
    temp2 = sumsq2/N - mean2*mean2
    if temp2 <= 0.:
        sig2 = 0.
    else:
        sig2 = np.sqrt(temp2)
    egFix = sumFixEst / N
    seMean.add(mean1, mean2, egFix, sig1/np.sqrt(N-1), sig2/np.sqrt(N-1))  # for each iTheta

    seRaw_data = data[2]
    seRaw.addEnsemble(seRaw_data)

    seNoFix_data = data[3]
    seNoFix.addEnsemble(seNoFix_data)

    dataf.close()
    #sys.exit()


print 'seMean %.6f, %.6f' % seMean.shear(), '+/- %.6f, %.6f' % seMean.sigmaE(shapeNoise=False)
print 'seRaw %.6f, %.6f' % seRaw.shear()
print 'seNoFix %.6f, %.6f' % seNoFix.shear()

"""
notetext = " "
P.figtext(0.52, 0.925, notetext, fontsize='xx-small', multialignment='left', weight='bold',)

annotatestr = S.Popen(['pwd'], stdout=S.PIPE).communicate()[0]
annotatestr = annotatestr + S.Popen(['date'], stdout=S.PIPE).communicate()[0]
annotatestr = annotatestr + sys.argv[0] + '::' + datafname
P.figtext(0.1, 0.93, annotatestr, fontsize='xx-small', multialignment='left',)

figfname = sys.argv[0].replace('.py', '.png')
P.savefig(figfname)
P.clf()
"""
