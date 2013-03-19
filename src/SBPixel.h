// $Id: SBPixel.h,v 1.5 2012/02/17 01:59:09 garyb Exp $
// SBProfile represented by interpolated pixel data table(s).
#ifndef SBPIXEL_H
#define SBPIXEL_H

#include "SBProfile.h"
#include "Interpolant.h"
#include "UseTMV.h"

namespace sbp {

  class SBPixel: public SBProfile {
  public:
    SBPixel(int Npix, double dx_, const Interpolant2d& i, int Nimages_=1);
    SBPixel(Image<> img, const Interpolant2d& i, double dx_=0., double padFactor=0.);
    SBPixel(const SBPixel& rhs);
    ~SBPixel();
    SBProfile* duplicate() const {return new SBPixel(*this);}

    // These are all the base class members that must be implemented:
    double xValue(Position<double> p) const;
    DComplex kValue(Position<double> p) const;
    // Notice that interpolant other than sinc may make max frequency higher than
    // the Nyquist frequency of the initial image
    double maxK() const {return xInterp->urange() * 2.*PI / dx;}
    // Require output FTs to be period on scale > original image extent + kernel footprint:
    double stepK() const {return 2.*PI / ( (Ninitial+2*xInterp->xrange())*dx);}
    bool isAxisymmetric() const {return false;}
    // This class will be set up so that both x and k domain values
    // are found by interpolation of a table:
    bool isAnalyticX() const {return true;}
    bool isAnalyticK() const {return true;}

    double centroidX() const;
    double centroidY() const;

    void setCentroid(Position<double> _p) {
      throw SBError("setCentroid not allowed for SBPixel");
    }

    double getFlux() const;
    void setFlux(double flux=1.);  // This will scale the weights vector

    /////////////////////
    // Methods peculiar to SBPixel

    // Set/get the value at one input pixel (without any weights or flux scaling)
    // Note that input data required to have -Npix/2 <= ix,iy < Npix/2
    // and 0 <= iz < Nimages.
    void setPixel(double value, int ix, int iy, int iz=0); //!! resets FTs
    double getPixel(int ix, int iy, int iz=0) const;

    // Set/get the weight vector applied to different planes
    void setWeights(const DVector& wts_); // ??? check dimensions first!
    DVector getWeights() const {return wts;}

    // Select interpolant used in x space
    void setXInterpolant(const Interpolant2d& interp_) {xInterp=&interp_; ready=false;}
    const Interpolant2d& getXInterpolant() const {return *xInterp;}

    // Select interpolant used in k space
    void setKInterpolant(const Interpolant2d& interp_) {kInterp=&interp_;}
    const Interpolant2d& getKInterpolant() const {return *kInterp;}

    // Return size of input data grid
    int getNin() const {return Ninitial;}
    // And size of DFT used to make k space:
    int getNft() const {return Nk;}

    // Overrides for better efficiency with separable kernels:
    virtual void fillKGrid(KTable &kt) const;
    virtual void fillXGrid(XTable &xt) const;
    virtual double fillXImage(Image<> I, double dx) const;

  private:
    void checkReady() const;	// Make sure all internal quantities are ok.
    int Ninitial;	// Size of input pixel grids
    double dx;		// and their pixel scales
    int Nk;		// Size of the padded grids and FT's
    double dk;		// step size in k tables
    int Nimages;	// Number of image planes to sum
    const Interpolant2d* xInterp;
    const Interpolant2d* kInterp;
    // set true if kTables, centroid/flux values,etc., are set for current x pixel values:
    mutable bool ready;	
    // Weights to use for summing each image plane:
    DVector wts;
    // Fluxes and x/y weighted fluxes for each image plane:
    mutable DVector fluxes;
    mutable DVector xFluxes;
    mutable DVector yFluxes;

    // The default k-space interpolant:
    static InterpolantXY defaultKInterpolant2d;

    // The input data arrays:
    vector<XTable*> vx;

    // The will need a bunch of mutable stuff for kTables and interpolations
    mutable vector<KTable*> vk;
    // Arrays summed with weights:
    mutable XTable* xsum;
    mutable KTable* ksum;
    mutable bool xsumValid;
    mutable bool ksumValid;
    void checkXsum() const;  // build xsum if it's not current
    void checkKsum() const;  // ...and ksums
  };
} // namespace sbp

#endif // SBPIXEL_H
