# A sample Makefile.  Edit to suit your platform.

# CXX can specify any paths to includes that are absolute and will be passed to subdirs
CXX = clang++
# OPTFLAGS will be exported for subdir makes
OPTFLAGS = -O3 -DASSERT -fPIC

# Note for compiling Boost.Python:
#
# Apple Python frameworks are under /System/Library/Frameworks/...
# User-installed Python frameworks go into /Library/Frameworks/...
#
# Since I've compiled my Boost.Python against /Library/Frameworks/EPD64.framework,
# and I also run the EPD64 version of Python,
# I need to make sure that all the includes and libs are from this framework.
# Otherwise, everything may compile, but when I load the fdnt module under python, I get the error:
# > TypeError: __init__() should return None, not 'NoneType'
# when instantiating C++ derived class.

# INCLUDES can be relative paths, and will not be exported to subdirectory makes.
INCLUDES = -I src -I src/utilities2 -I src/images -I src/astrometry2 -I pysrc \
           -I /Library/Frameworks/EPD64.framework/Versions/7.2/include/python2.7 \
           -I /Library/Frameworks/EPD64.framework/Versions/7.2/lib/python2.7/site-packages/numpy/core/include \
           -I /opt/local/include \
           -I /Users/reiko/2code/cfitsio/include

CXXFLAGS = $(OPTFLAGS) $(INCLUDES)
SRC = $(shell ls *.cpp)

SUBDIRS = src src/utilities2 src/images src/astrometry2

TMV_LINK := $(shell cat /usr/local/share/tmv/tmv-link)

LIBS = -lm \
	-lboost_python \
	-L/Library/Frameworks/EPD64.framework/Versions/7.2/lib -lpython2.7 \
	-L/opt/local/lib -lfftw3 \
	-L/Users/reiko/2code/cfitsio/lib -lcfitsio \
	-ltmv_symband $(TMV_LINK) \
	-L/Users/reiko/2code/CCfits/.libs -lCCfits \
	-L/usr/local/lib -ltmv_symband -ltmv \
	-lpthread

SUBOBJ = src/utilities2/BinomFact.o \
	src/utilities2/Interpolant.o \
	src/utilities2/StringStuff.o \
	src/utilities2/fft.o \
	src/utilities2/GTable.o \
	src/utilities2/Pset.o \
	src/utilities2/Poly2d.o \
	src/images/FITS.o \
	src/images/Image.o \
	src/images/FITSImage.o \
	src/images/HeaderFromStream.o \
	src/astrometry2/PixelMap.o \
	src/astrometry2/Astrometry.o \
	src/astrometry2/PolyMap.o \
	src/astrometry2/PixelMapCollection.o

OBJ = src/Laguerre.o \
	src/Shear.o \
	src/LTransforms.o \
	src/GLSimple.o \
	src/SBProfile.o \
	src/SBPixel.o \
	src/SCAMPMap.o \
	src/PSFEx.o \
	src/FDNT.o \
	src/ExposureGroup.o \
	src/Lance.o $(SUBOBJ)

all: depend subs

_fdnt: module.os FDNTImage.os RunFDNT.os RunFDNT.o Bounds.os
	$(CXX) -bundle pysrc/.obj/module.os pysrc/.obj/FDNTImage.os pysrc/.obj/RunFDNT.os pysrc/.obj/Bounds.os src/.obj/RunFDNT.o $(OBJ) \
	$(LIBS) -o fdnt/$@.so

module.os: pysrc/module.cpp
	$(CXX) $(CXXFLAGS) $^ -c -o pysrc/.obj/$@

FDNTImage.os: pysrc/FDNTImage.cpp
	$(CXX) $(CXXFLAGS) $^ -c -o pysrc/.obj/$@

RunFDNT.os: pysrc/RunFDNT.cpp
	$(CXX) $(CXXFLAGS) $^ -c -o pysrc/.obj/$@

Bounds.os: pysrc/Bounds.cpp
	$(CXX) $(CXXFLAGS) $^ -c -o pysrc/.obj/$@

RunFDNT.o: src/RunFDNT.cpp
	$(CXX) $(CXXFLAGS) $^ -c -o src/.obj/$@

FDNTColorGrad: src/FDNTColorGrad.o src/ColorGrad.o src/EnclosedFluxRadius.o SBParse.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o bin/$@
FDNTPSFEx: src/FDNTPSFEx.o src/EnclosedFluxRadius.o src/SBParse.o src/images/HeaderFromStream.o \
	$(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o bin/$@
FDNTGrid: src/FDNTGrid.o src/EnclosedFluxRadius.o src/SBParse.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o bin/$@
FDNTRing_spl: src/FDNTRing_spl.o src/EnclosedFluxRadius.o src/SBParse.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o bin/$@
FDNTRing_omp: src/FDNTRing_omp.o src/EnclosedFluxRadius.o src/SBParse.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o bin/$@
FDNTRing: src/FDNTRing.o src/EnclosedFluxRadius.o src/SBParse.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o bin/$@
FDNTEllipse: src/FDNTEllipse.o src/EnclosedFluxRadius.o src/SBParse.o src/Lance.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o bin/$@
FDNTGreat: src/FDNTGreat.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o bin/$@

###############################################################
## Standard stuff:
###############################################################

export OPTFLAGS
export CXX

subs:
	for dir in $(SUBDIRS); do (cd $$dir && $(MAKE)); done

depend:
	for dir in $(SUBDIRS); do (cd $$dir && $(MAKE) depend); done
# there is no .cpp file in the head node
        #$(CXX) $(CXXFLAGS) -MM $(SRC) > .$@

clean:
	for dir in $(SUBDIRS); do (cd $$dir && $(MAKE) clean); done
	rm -f src/*.o src/*~ src/*.dvi src/*.aux core .depend

ifeq (.depend, $(wildcard .depend))
include .depend
endif

.PHONY: all install dist depend clean 