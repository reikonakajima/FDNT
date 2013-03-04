#ifndef PSFEX_H
#define PSFEX_H

#include "SBPixel.h"

namespace sbp {

  class PSFExError: public SBError {
  public:
    PSFExError(const string &m=""): SBError(m) {}
  };

  class PSFExModel
  {
  public:
    enum PSFExFormat {
      TEXT,
      FITS
    };
     
    PSFExModel(const char *filename, PSFExFormat format=TEXT);
    PSFExModel(const PSFExModel& rhs);
    ~PSFExModel();

    const SBPixel* sb() const {return sbp;}
    void fieldPosition(double x, double y);
    // Note that setting flux will be un-done whenever new field 
    // position is set.
    void setFlux(double flux) {sbp->setFlux(flux);}

    // Set this to use cubic interpolant in k space instead of default quintic if you're
    // worried about speed.
    void useFasterKInterpolant(bool useCubic=true) {
      sbp->setKInterpolant(useCubic ? cubic_2d : quintic_2d);
    }

    // return size of model pixels    
    double getDx() const {return dx;}
    // Get order of interpolating polynomial
    int getOrder() const {return order;}
    // A number provided by PSFEx:
    double getFWHM() const {return fwhm;}

  private:
    int polnaxis;  // Note we are only accepting 2-axis polynomials
    // Information on each axis of the polynomial:
    vector<int> polgrp;  // Only working with polgrp=1 in this code
    vector<double> polzero;
    vector<double> polscal;
    vector<string> polname;

    vector<int> poldeg; // degree of each polynomial group

    double fwhm;
    double dx;  // "sampling" in PSFEx language

    int Nx;	// Larger of x and y axis lengths
    int order;
    SBPixel* sbp;

    // order in field x and y of polynomial term multiplying each image plane:
    vector<int> xOrder;
    vector<int> yOrder;

    static InterpolantXY lanczos3_2d;
    static InterpolantXY cubic_2d;
    static InterpolantXY quintic_2d;
  };

} // namespace sbp
#endif // PSFEX_H
