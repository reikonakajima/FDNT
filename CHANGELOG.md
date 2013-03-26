
Initial Repository: modifications from original FDNT code:

#1
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


#2
==
Owner:  Reiko
Topic:  Various modifications
	- Modify FDNTRing as FDNTRing_omp to allow multi-processing using OpenMP
	- Modify FDNTRing as FDNTRing_spl to run piecewise on multiple processes
	- Made no size marginalization default behaviour of FDNT.  This allows the code to 
	  run ~10x faster, for high S/N cases.
        - Tested: Smaller weight does not help with reducing scatter for undersampled Airy PSFs.

#3
==
Owner:  Ole
Topic:  Debug FDNTGreat gridding bug: fix on the grid assignment
