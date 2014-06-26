Installation Instructions
=========================

System requirements
-------------------

SCons must be installed.
A C++ compiler, such as g++, is required.


Required software
-----------------

Install the following libraries

  - CCfits
  - cfitsio
  - tmv v0.72  or later   from Mike Jarvis
  - fftw v3.0  or later
  - Python v2.x (tested with v.2.7)
  - boost_python, compiled against the python in use



SCons options (relevant after issue #6 mergers with master)
-----------------------------------------------------------

The SConstruct/SConscript file has been heavily copied from that of the GalSim version.

To see all the SCons options, type `scons --help`.

The most relevant options to set are:
PREFIX               Your final installation directory
PYPREFIX             Your final python module installation directory
TMV_DIR              Your tmv installation directory
BOOST_DIR            Your Boost installation directory
EXTRA_INCLUDE_PATH   Where you can find other headers (such as cfitsio.h)
EXTRA_LIB_PATH       Where you can find other libs (such as libCCfits.h)




Edit the main Makefile (soon to be obsolete)
--------------------------------------------

There is a Makefile template file named `src/sample.Makefile`.  Make a copy of 
this file and name it `src/Makefile`.  Ensure that the includes and lib files 
are properly reached in the Makefile options.

The `src` contains subdirectories `astrometry2`, `images`, and `utilities2`, 
each of which contain its own Makefile.  These need not be edited.

Note that the file `.gitignore` contains the line `Makefile*`.  This means git 
will ignore any files starting with letters `Makefile`---i.e., it will leave
your system-specific `Makefile` alone.
