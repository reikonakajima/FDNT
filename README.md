FDNT
====

Fourier Domain Null Test (FDNT) galaxy shape measurement code.


About the code
--------------

* The main code is written in C++.  Some scripts are also written in Python.

* The code is not doing exactly what's in the FDNT paper (Bernstein 2010).
  In particular:
  - the criterion for setting the weight size is changed.
  - the baseline method is to map the likelihood space in the (eta1,eta2)
    plane and find the mean of the distribution, rather than searching for
    maximum likelihood.  
  - this code has the remnants of various experiments in it. 

* The key calculations for the FDNT method are in the tests() routine in
  ExposureGroup.cpp, so that might be a place to start trying to
  comprehend. 

* For installation instructions, consult the file INSTALL.md.

* The various update history to the code should be described in the file
  CHANGELOG.md.

* Please note the condition of usage of this code, in the WARNING.md file.



How to communicate with the FDNT developers
-------------------------------------------

If you have any comments, questions, or suggestions, please open up an Issue on
our GitHub repository:

https://github.com/reikonakajima/FDNT/issues?state=open



Installation
------------

For installation instructions, please see the file `INSTALL.md` in the main
repository directory.



Repository directory structure
------------------------------

The repository has a number of subdirectories. Below is a guide to their
contents:

* bin/ :      executables (after the compilation procedure is done).
* devel/ :    an assortment of developer tools.
* doc/ :      documentation on differnt portions of the code.
* src/ :      the source C++ code of FDNT.
* pysrc/ :    Python scripts.
* pipeline/ : a pipeline setup using FDNTGreat.

