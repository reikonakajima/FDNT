#   $Id: Makefile,v 1.6 2012/02/02 23:34:25 garyb Exp $ 
SUBDIRS = 
UTIL = ../utilities2
INCLUDES = -I $(UTIL)
SRC = $(shell ls *.cpp)
OPTFLAGS = -O3
CPPFLAGS = $(OPTFLAGS) $(INCLUDES)
OBJ = Astrometry.o Ephemeris.o PixelMap.o PolyMap.o PixelMapCollection.o

all: depend $(OBJ)

SUBOBJ = $(UTIL)/StringStuff.o 

depend:
	for dir in $(SUBDIRS); do (cd $$dir && $(MAKE) depend); done
	$(CXX) $(CPPFLAGS) -MM $(SRC) > .$@

clean:
	for dir in $(SUBDIRS); do (cd $$dir && $(MAKE) clean); done
	rm -f *.o *~ core .depend

ifeq (.depend, $(wildcard .depend))
include .depend
endif

export
# don't export INCLUDES??

.PHONY: all install dist depend clean 
