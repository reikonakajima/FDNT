// $Id: Lance.cpp,v 1.3 2011/04/13 03:43:36 garyb Exp $
#include "Lance.h"
#include <cmath>

// Implement hexagonal-grid points:  here is
// xy notation:

//    25      35     45      55     65

// 24     34      44     54     64

//    23      33     43      53     63

// 22     32      42     52     62

//    21      31     41      51     61

// 20     30      40     50     60

Sample::Sample(double x0, double y0, double step_, int level_):
  xOrigin(x0), yOrigin(y0), step(step_), level(level_), 
  ix(0), iy(0), lnProb(-999.), sigma(-999.), covE(2) {}

bool
Sample::operator<(const Sample& rhs) const {
  return iy < rhs.iy || ( iy==rhs.iy && ix<rhs.ix) ;
}

void
Sample::getXY(double& x, double& y) const {
  const double yFactor = std::sqrt(3.)/2.;
  y = yOrigin + iy*yFactor*step;
  // For hex, odd-numbered rows are shifted right by 1/2 step
  x = xOrigin + (ix + (iy%2==0 ? 0. : 0.5))*step;
}


list<Sample>
Sample::neighbors(int outLevel) const {
  list<Sample> ls;
  // Return empty list if attempting below base grid resolution
  if (outLevel<0) return ls;

  Sample out=*this;
  out.level = outLevel;

  // Distance to neighbors is 2^(outLevel);
  int istep=1;
  for (int i=0; i<outLevel; i++) istep*=2;
  // Add point to the right:
  out.ix = ix + istep;
  ls.push_back(out);
  // Add point to the left:
  out.ix = ix - istep;
  ls.push_back(out);
  // Above left:
  out.iy = iy + istep;
  if (level==1)
    out.ix = ix - (iy%2==0 ? 1 : 0);
  else
    out.ix = ix - istep/2;
  ls.push_back(out);
  // Below left:
  out.iy = iy - istep;
  ls.push_back(out);
  // Below right
  out.ix += istep;
  ls.push_back(out);
  // Above right
  out.iy = iy + istep;
  ls.push_back(out);

  return ls;
}

list<Sample>
Sample::ringAround(int radius) const {
  list<Sample> ls;
  // Nothing for null radius
  if (radius<=0) return ls;

  Sample out=*this;

  // Grid step is 2^level:
  int istep=1;
  for (int i=0; i<level; i++) istep *= 2;

  // Start to the right of this point:
  out.ix = ix + radius * istep;
  ls.push_back(out);
  // Move up & to left radius times;
  for (int i=0; i<radius; i++) {
    if (level>0) out.ix-=istep/2;
    else out.ix -= (out.iy%2==0 ? 1 : 0);
    out.iy += istep;
    ls.push_back(out);
  }
  // Move left:
  for (int i=0; i<radius; i++) {
    out.ix -= istep;
    ls.push_back(out);
  }
  // Move down & to left
  for (int i=0; i<radius; i++) {
    if (level>0) out.ix-=istep/2;
    else out.ix -= (out.iy%2==0 ? 1 : 0);
    out.iy -= istep;
    ls.push_back(out);
  }
  // Move down & to right
  for (int i=0; i<radius; i++) {
    if (level>0) out.ix+=istep/2;
    else out.ix += (out.iy%2==0 ? 0 : 1);
    out.iy -= istep;
    ls.push_back(out);
  }
  // Move right
  for (int i=0; i<radius; i++) {
    out.ix += istep;
    ls.push_back(out);
  }
  // Move up & to right - already have last point
  for (int i=0; i<radius-1; i++) {
    if (level>0) out.ix+=istep/2;
    else out.ix += (out.iy%2==0 ? 0 : 1);
    out.iy += istep;
    ls.push_back(out);
  }
  return ls;
}

    
