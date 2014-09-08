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
CCFITS_DIR           Your CCfits installation directory
EXTRA_INCLUDE_PATH   Where you can find other headers (such as cfitsio.h)
EXTRA_LIB_PATH       Where you can find other libs


The compile and install commands are:

   >>> scons --PREFIX=/your/chosen/directory
   >>> scons install

or

   >>> scons
   >>> sudo scons install

To remove installation:

   >>> scons -c
   >>> scons -c install  (or sudo scons -c install)

Sometimes, the a new compilation will fail.  Try

   >>> /bin/rm -r .scon*

before re-running scons compilation.
