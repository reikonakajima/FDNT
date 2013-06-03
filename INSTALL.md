Installation Instructions
=========================

System requirements
-------------------

A C++ compiler, such as gcc, is required.


Required software
-----------------

Install the following libraries

  - CCfits
  - cfitsio
  - tmv v0.70  or later   from Mike Jarvis
  - fftw v3.0  or later


Edit the main Makefile
----------------------

There is a Makefile template file named `src/sample.Makefile`.  Make a copy of 
this file and name it `src/Makefile`.  Ensure that the includes and lib files 
are properly reached in the Makefile options.

The `src` contains subdirectories `astrometry2`, `images`, and `utilities2`, 
each of which contain its own Makefile.  These need not be edited.

Note that the file `.gitignore` contains the line `Makefile*`.  This means git 
will ignore any files starting with letters `Makefile`---i.e., it will leave
your system-specific `Makefile` alone.