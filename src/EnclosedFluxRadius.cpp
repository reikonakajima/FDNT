// $Id: EnclosedFluxRadius.cpp,v 1.2 2011/07/20 15:44:35 dgru Exp $
// Report the radius enclosing desired flux fraction
// for an Image or an SBProfile.

#include <vector>
#include <cstdlib>
#include <string>
#include <cmath>
#include <stdexcept>

#include "SBProfile.h"
#include "Image.h"

using namespace std;

class Pair {
public:
  Pair(double x_, double y_): x(x_), y(y_) {}
  double x, y;
  bool operator<(const Pair& rhs) const {return x<rhs.x;}
};

double EnclosedFluxRadius(Image<> img, double enclosedFraction=0.5) {
  Bounds<int> b=img.getBounds();

  double xsum=0., ysum=0., flux=0.;
  for (int j=b.getYMin(); j<=b.getYMax(); j++)
    for (int i=b.getXMin(); i<=b.getXMax(); i++) {
      xsum += i*img(i,j);
      ysum += j*img(i,j);
      flux += img(i,j);
    }

  double xc = xsum/flux;
  double yc = ysum/flux;

  list<Pair> lp;
  for (int j=b.getYMin(); j<=b.getYMax(); j++)
    for (int i=b.getXMin(); i<=b.getXMax(); i++) {
      double rad=hypot(i-xc, j-yc);
      double val=img(i,j);
      lp.push_back(Pair(rad,val));
    }
  lp.sort();
  double fsum=0.;
  double eer=0.;
  for (list<Pair>::iterator i=lp.begin();
       i != lp.end();
       ++i) {
    double val = i->y;
    fsum += val;
    if (fsum >= enclosedFraction*flux & eer==0.) {
	// Fit included flux to a line in this vicinity
	double sxx, sxy, syy, sx, sy;
	double ff=fsum;
	double y=ff-enclosedFraction*flux;
	double x0 = i->x;
	double x=0.;
	sxx = x*x; sxy = x*y; syy = y*y; 
	sx = x; sy=y;
	int n=1;
	list<Pair>::iterator j=i;
	// Go down:
	while ((n<=3 || j->x > 0.9*i->x) && j!=lp.begin()) {
	  ff -= j->y;
	  j--;
	  y=ff-enclosedFraction*flux;
	  x=j->x-x0;
	  sxx += x*x; sxy += x*y; syy += y*y; 
	  sx += x; sy+=y;
	  n++;
	}
	// Go up:
	j = i;
	int nup=1;
	ff = fsum;
	for (j++; j!=lp.end() && (nup<=2 || j->x<1.1*i->x); j++, nup++) {
	  ff += j->y;
	  y=ff-enclosedFraction*flux;
	  x = j->x-x0;
	  sxx += x*x; sxy += x*y; syy += y*y; 
	  sx += x; sy+=y;
	  n++;
	}
	// intercept & slope formula:
	double m = (n*sxy - sx*sy)/(n*sxx-sx*sx);
	double b = (sxx*sy - sx*sxy)/(n*sxx-sx*sx);
	return -b/m + x0;
      }
    }
  // Should not get here
  throw ImageError("EnclosedFluxRadius not reached");
}

// Find EEradius for an SBProfile
double EnclosedFluxRadius(const sbp::SBProfile& sbp, double enclosedFraction=0.5) {
  // First let it choose its own dx for drawing
  double dx;
  double eer;
  {
    Image<> img = sbp.draw();
    img.getHdrValue("DX",dx);
    eer = EnclosedFluxRadius(img, enclosedFraction);
  }
  // eer is now in pixels.  Want a minimum size for accuracy:
  const double minimumEEpix = 6.;
  const double targetEEpix = 10.;
  if (eer < minimumEEpix) {
    dx = eer*dx / targetEEpix;
    Image<> img2 = sbp.draw(dx);
    img2.getHdrValue("DX",dx);
    eer = EnclosedFluxRadius(img2, enclosedFraction);
    if (eer < minimumEEpix)
      throw ImageError("EnclosedFluxRadius could not draw SBProfile large enough");
  }
  return eer*dx;
}
