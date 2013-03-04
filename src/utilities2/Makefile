#   $Id: Makefile,v 1.3 2012/01/01 18:56:10 garyb Exp $ 
SUBDIRS = 
OPTFLAGS = -O3
CXXFLAGS = $(OPTFLAGS) -I /usr/local/include -I /sw/include
SRC = $(shell ls *.cpp)
OBJ = BinomFact.o fft.o StringStuff.o Poisson.o Table.o Pset.o odeint.o Poly2d.o Interpolant.o

all: depend $(OBJ)

depend:
	for dir in $(SUBDIRS); do (cd $$dir && $(MAKE) depend); done
	$(CXX) $(CXXFLAGS) -MM $(SRC) > .$@

clean:
	for dir in $(SUBDIRS); do (cd $$dir && $(MAKE) clean); done
	rm -f *.o *~ core .depend

ifeq (.depend, $(wildcard .depend))
include .depend
endif

export

.PHONY: all install dist depend clean 
