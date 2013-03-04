
This file keeps track of the various branches in the 'fdnt1' repository.


#5_Test_Airy_Weight
===================
Owner:  Reiko
Topic:  Airy PSF with DeVauc is unstable, see if smaller weight would help
Start:  Mar 4, 2013
Status: Coding
Merged: -


#4_Verify_FDNTPSFEx
===================
Owner:  ? (Ole or Reiko)
Topic:  Verify that FDNTPSFEx (merged from DG into master) works
Start:  Mar 4, 2013
Status: Hasn't started yet
Merged: -


#3_FDNTRing_no_size_marginalization
===================================
Owner:  Reiko Nakajima
Topic:  Generate and test code without size marginalization
Start:  Feb 27, 2013
Status: Behavior confirmed on high S/N cases
Merged: Mar 4, 2013


DG
==
Owner:  Daniel Gruen
Topic:  Various modifications
        - assume that the initial guess of ellipticity is actually a
       	  'good guess' of pre-seeing ellipticity (i.e., not strongly
	  influenced by PSF size and shape as the SExtractor ellipse
	  I had used previously).
        - In the clusterSTEP simulations, I ended up using a KSB shape
       	  estimate as the starting point for sampling probability space.
        - some of them are quite quick&dirty fixes of code which caused
          crashes in rare occasions on real data.
        - [The repository does not yet include the recent (or, for that matter,
       	  any) version of an FDNT code that can read PSFEx models of the PSF.]
Start:  Aug ??, 2012 or prior
Status: Behavior confirmed on high S/N cases
Merged: Mar 4, 2013


#2_FDNTRing_imperial
====================
Owner:  Reiko
Topic:  Modify FDNTRing to run piecewise on multiple processes (on Imperial)
Start:  Feb 19, 2013
Status: Verified at high S/N
Merged: Feb 21, 2013


#1
==
Owner:  Ole
Topic:  Debug FDNTGreat gridding bug
Start:  Feb 18, 2013
- 20130218: Initial. Clone of master, plus pipeline prototype from RN v2 repo
Status: Verifying on GREAT08-like images, waiting for another party (reiko) to finish check
Merged: -


#1_FDNTRing_omp
===============
Owner:  Reiko
Topic:  Allow multi-processing using OMP of FDNTRing tests
Start:  Feb ??, 2013
Status: Verified at high S/N
Merged: Feb 19, 2013


lnsigstep
=========
Owner:  Probably Reiko
Topic:  Modify step size of ln(sigma), the size of the galaxy, for marginalizing
Start:  Aug ??, 2012 (probably)
Status: (probably taken over by master)
Merged: (no, leave it alone)
