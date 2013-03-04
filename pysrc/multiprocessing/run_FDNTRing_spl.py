#!/usr/bin/env python

import subprocess as S
import numpy as np
import pylab as P
import sys
sys.path.append('/Users/reiko/2code/python/mylib')

"""
run_FDNTRing_spl.py
  allow distribution of FDNTRing tests over many cores.

"""

command = '../src/FDNTRing_spl'
nTheta = 32
runCondor = True
if runCondor:
    # generate condor file
    condorfilename = 'test.condor'
    condorf = open(condorfilename, 'w')


shcommand = '/bin/rm *.out *.par *.log *.err'
S.call(shcommand, shell=True)

first_line = True
for iTheta in range(nTheta):

    tag = 'test' + str(iTheta)
    print tag

    arguments = """sbPSF = (airy 0.812 0.33 S 0.1 0.1) * box
sbGalaxy = exp 1.4 S 0.3 0.1
nTheta = %d
iTheta = %d
nDither = 3
nNoise = 10
sig = 10000""" % (nTheta, iTheta,)

    argfile = tag + '.par'
    argf = open(argfile, 'w')
    print >> argf, arguments
    argf.close()

    if not runCondor:
        outfile = tag + '.out'
        shcommand = '%s < %s > %s' % (command, argfile, outfile)
        S.Popen(shcommand, shell=True)   # let them all run in parallel !!  (forks)
    
    else:
        if first_line:
            header = """Universe    = vanilla
Executable  = /home/reiko/2code/fdnt1/src/FDNTRing_spl
Arguments   = /home/reiko/2code/fdnt1/pysrc/%s
Log         = test.log
Output      = %s.out
Error       = %s.err
Queue""" % (argfile, tag, tag,)
            print >> condorf, header
            first_line = False
        else:
            repeater = """
Arguments   = /home/reiko/2code/fdnt1/pysrc/%s
Output      = %s.out
Error       = %s.err
Queue""" % (argfile, tag, tag,)
            print >> condorf, repeater

condorf.close()



# sbPSF = airy 0.3155 0.33 * box

# sbPSF = (airy 0.812 0.33 S 0.1 0.1) * box  # AiryR09
