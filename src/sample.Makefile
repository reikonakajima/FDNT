# $Id: Makefile,v 1.20 2012/01/02 19:25:58 garyb Exp $
# make the test programs

# CXX can specify any paths to includes that are absolute and will be passed to subdirs
CXX = g++ -fopenmp
# OPTFLAGS will be exported for subdir makes
OPTFLAGS = -O3 -DASSERT

# INCLUDES can be relative paths, and will not be exported to subdirectory makes.
INCLUDES = -I utilities2 -I images -I astrometry2 -I /euclid/data02/euclid01/local/include

CXXFLAGS = $(OPTFLAGS) $(INCLUDES)
SRC = $(shell ls *.cpp)

SUBDIRS = utilities2 images astrometry2

TMV_LINK := $(shell cat /euclid/data02/euclid01/local/share/tmv/tmv-link)

LIBS = -lm -L/opt/local/lib -lfftw3 -lcfitsio -ltmv_symband $(TMV_LINK) -L/usr/lib/x86_64-linux-gnu/ -lCCfits -L/euclid/data02/euclid01/local/lib

SUBOBJ = utilities2/BinomFact.o images/FITS.o utilities2/Interpolant.o utilities2/BinomFact.o \
	utilities2/StringStuff.o images/Image.o images/FITSImage.o images/HeaderFromStream.o \
	utilities2/fft.o utilities2/GTable.o utilities2/Pset.o utilities2/Poly2d.o \
	astrometry2/PixelMap.o astrometry2/Astrometry.o astrometry2/PolyMap.o \
	astrometry2/PixelMapCollection.o

OBJ = Laguerre.o Shear.o LTransforms.o GLSimple.o SBProfile.o SBPixel.o SCAMPMap.o PSFEx.o \
	FDNT.o ExposureGroup.o Lance.o $(SUBOBJ)

all: depend subs

PSFEx_KiDS_astrom_radec_only: PSFEx_KiDS_astrom_radec_only.o EnclosedFluxRadius.o SBParse.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o ../bin/$@
PSFEx_KiDS_astrom: PSFEx_KiDS_astrom.o EnclosedFluxRadius.o SBParse.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o ../bin/$@
PSF_KiDS: PSF_KiDS.o EnclosedFluxRadius.o SBParse.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o ../bin/$@
PSFEx_KiDS: PSFEx_KiDS.o EnclosedFluxRadius.o SBParse.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o ../bin/$@
FDNTPSFEx_GREAT3_varPSF: FDNTPSFEx_GREAT3_varPSF.o EnclosedFluxRadius.o SBParse.o images/HeaderFromStream.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o ../bin/$@
FDNTPSFEx_GREAT3_debug: FDNTPSFEx_GREAT3_debug.o EnclosedFluxRadius.o SBParse.o images/HeaderFromStream.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o ../bin/$@
spacePSF_debug: spacePSF_debug.o EnclosedFluxRadius.o SBParse.o images/HeaderFromStream.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o ../bin/$@
FDNTPSFEx_GREAT3: FDNTPSFEx_GREAT3.o EnclosedFluxRadius.o SBParse.o images/HeaderFromStream.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o ../bin/$@
FDNTRing_spl_true_ic: FDNTRing_spl_true_ic.o EnclosedFluxRadius.o SBParse.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o ../bin/$@
FDNTColorGrad: FDNTColorGrad.o ColorGrad.o EnclosedFluxRadius.o SBParse.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o ../bin/$@
FDNTPSFEx: FDNTPSFEx.o EnclosedFluxRadius.o SBParse.o images/HeaderFromStream.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o ../bin/$@
FDNTGrid: FDNTGrid.o EnclosedFluxRadius.o SBParse.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o ../bin/$@
FDNTRing_spl: FDNTRing_spl.o EnclosedFluxRadius.o SBParse.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o ../bin/$@
FDNTRing_omp: FDNTRing_omp.o EnclosedFluxRadius.o SBParse.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o ../bin/$@
FDNTRing: FDNTRing.o EnclosedFluxRadius.o SBParse.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o ../bin/$@
FDNTEllipse: FDNTEllipse.o EnclosedFluxRadius.o SBParse.o Lance.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o ../bin/$@
FDNTGreat: FDNTGreat.o $(OBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o ../bin/$@

###############################################################
## Standard stuff:
###############################################################

export OPTFLAGS
export CXX

subs:
	for dir in $(SUBDIRS); do (cd $$dir && $(MAKE)); done

depend:
	for dir in $(SUBDIRS); do (cd $$dir && $(MAKE) depend); done
	$(CXX) $(CXXFLAGS) -MM $(SRC) > .$@

clean:
	for dir in $(SUBDIRS); do (cd $$dir && $(MAKE) clean); done
	rm -f *.o *~ *.dvi *.aux core .depend

ifeq (.depend, $(wildcard .depend))
include .depend
endif

.PHONY: all install dist depend clean 
