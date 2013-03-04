// $Id: Lance.h,v 1.3 2011/04/13 03:43:36 garyb Exp $
// Classes to implement variant of Lance Miller's e-plane search algorithm

#include <list>
#include "UseTMV.h"
#include "TMV_SymBand.h"
using std::list;

#ifndef LANCE_H
#define LANCE_H
class Sample {
public:
  // Constructor builds node at origin
  Sample(double x0, double y0, double step, int level);
  ~Sample() {}
  bool operator<(const Sample& rhs) const;

  // Return all points at current level lying on
  // a hexagonal "ring" of size radius * (2^level)
  list<Sample> ringAround(int radius) const;

  // Return the six hex neighbors centered on current
  // point at resolution outLevel 
  list<Sample> neighbors(int outLevel) const;

  // No arguments returns neighbors at current level
  list<Sample> neighbors() const {return neighbors(level);}
  // spawn() returns neighbors at daughter level (densification)
  list<Sample> spawn() const {return neighbors(level-1);}
  void getXY(double& x, double& y) const;

  int level;
  double lnProb;

  // More information about likelihood calc here:
  double sigma;
  double fractionTooBig;
  double fractionTooSmall;
  bool centroidFlag;
  tmv::SymMatrix<double> covE;
private:
  int ix;
  int iy;
  double xOrigin;
  double yOrigin;
  double step;
};

#endif   // LANCE_H
