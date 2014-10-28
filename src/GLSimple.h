// $Id: GLSimple.h,v 1.10 2011/11/04 18:20:29 dgru Exp $
//
// A class definition for fitting - fitting on multiple images is possible.
//  
#ifndef GLSIMPLE_H
#define GLSIMPLE_H

#include <cmath>
#include <list>

//#define DEBUGLOGGING
#include "Std.h"
#include "UseTMV.h"
#include "TMV_SymBand.h"
#include "Shear.h"
#include "Image.h"
#include "Laguerre.h"
#include "FitFlags.h"

namespace laguerre {

  class GLError: public std::runtime_error {
  public:
    GLError(const string &m=""): std::runtime_error("GL Error: " +m) {}
  };

  // A little class used to map world coords back into pixels,
  // roughly, when getting new masks.
  class CrudeMap {
  private:
    // Map is defined as w-w0 = dWdP * (p-p0);
    DMatrix dWdP;	
    DMatrix dPdW;	
    DVector w0;
    DVector p0;
  public:
    void pixToWorld(double px, double py, double& wx, double& wy);
    void worldToPix(double wx, double wy, double& px, double& py);
    Ellipse worldToPix(const Ellipse& ew);
    // Initialize by getting local map from point in an image close to
    // the specified world coord system.  Use near edge if off image.
    CrudeMap(const Image<double> xw, const Image<double> yw,
	     double xw0, double yw0);
    DMatrix getdPdW() const {return dPdW;}
  };

  // Structure to hold the information needed to fit a chosen exposure
  // Note all elements are pointers, fundamental types, or have pointer semantics,
  // so this structure can be passed around freely and all creation/destruction may help.
  template <class T=float>
  class FitExposure {
  public:
    const Image<T> sci;
    const Image<T> wt;   // inverse variance
    Image<double> xworld;	// world coordinates at each input pixel
    Image<double> yworld;
    //conversion factor from sci units to surface brightness:
    const double sbScaleFactor;
    // Sky level in units of sci image:
    const double sky;
    // Identifies PSF information associated with this exposure - not
    // used in the native GLSimple routines;
    // it is used for mapping to the PSF list in FDNT and should be in [0,npsfs-1] for this purpose
    const int psfCode;	

    // Simple constructor with just data & error, pixel units = world
    FitExposure(const Image<T> sci_, const Image<T> wt_, 
		int psf_=-1, double sky_=0., double scale=1.):
      sci(sci_), wt(wt_), xworld(sci_.getBounds()), yworld(sci_.getBounds()),
      sbScaleFactor(scale), sky(sky_), psfCode(psf_)
    {
      for  (int iy=yworld.YMin(); iy<=yworld.YMax(); iy++)
	for (int ix=xworld.XMin(); ix<=xworld.XMax(); ix++) {
	  xworld(ix,iy)=ix; yworld(ix,iy)=iy;
	}
    }

    // Provide precomputed world coord images
    FitExposure(const Image<T> sci_, const Image<T> wt_, 
		Image<double> xworld_, Image<double> yworld_,
		const int psf_=-1, double sky_=0., double scale=1.):
      sci(sci_), wt(wt_), xworld(xworld_), yworld(yworld_),
      sbScaleFactor(scale), sky(sky_), psfCode(psf_) {}
      
    vector< Position<int> > badPixels(double &meanweight)
    {
      meanweight=0;
      vector< Position<int> > n;
      for  (int iy=wt.YMin(); iy<=wt.YMax(); iy++)
	for (int ix=wt.XMin(); ix<=wt.XMax(); ix++)
	{
	  if(wt(ix,iy)<=0.) n.push_back(Position<int>(ix,iy));
	  else meanweight += wt(ix,iy);
	}
      meanweight /= double((wt.YMax()-wt.YMin()+1)*(wt.XMax()-wt.XMin()+1)-n.size());
      return n;
    }
  };
  

  template <class T=float>
  class GLSimple {
  public:
    // Constructor for plain-old fitting, single image, pixel coords
    GLSimple(const Image<T> sci, const Image<T> wt,
	  Ellipse e_, int startingOrder=8, double sky_=0.);
    // Constructor with multiple exposures:
    GLSimple(list<FitExposure<T> > list_,
	     Ellipse e_, int startingOrder=8);
    // Constructor for one FitExposure:
    GLSimple(FitExposure<T>& fe,
	     Ellipse e_, int startingOrder=8);
    ~GLSimple();


    ///////////////////////////////////////////////
    // The routines of interest:
    // Return chisq-minimizing b vector at
    // currently specified order and basis.
    // Masks are redrawn if basis has changed since last fit.
    LVector linearSolution();

    // An iterative fit to satisfy the centroid/shear/size tests.
    // Uses currently specified order.
    // Return true on success
    bool solve();

    // Force re-masking of data (otherwise is done only during solve()):
    // Return is false if center is off all images or too large.
    bool reMask();

    ///////////////////////////////////////////////
    // Routines to control the behavior of the next fit:
    void setBasis(Ellipse e);
    const Ellipse& getBasis() const {return basis;}
    void setOrder(int ord_);
    int    getOrder() const {return order;}

    // Mask parameter changes do *not* force remasking: 
    // just affect future mask choices.
    void setMaximumMaskMu(double max)    {maxMaskMu=max;}
    double getMaximumMaskMu() {return maxMaskMu;}
    //Set to <=0. for automatic order-based choice
    void setMaskSigma(double sig) {maskSigma=sig; invalidateData();}

    bool getCentering() const {return centering;}
    void setCentering(bool f);
    bool getDilating() const {return dilating;}
    void setDilating(bool f);
    bool getShearing() const {return shearing;}
    void setShearing(bool f);

    void setTolerance(double tol) {tolerance=tol; isSolved=false;}
    void setMaxIterations(int max)   {maxIterations=max;}

    ///////////////////////////////////////////////
    // Return results of fits:

    int  getFlags() const {return flags;}
    void  clearFlags() {flags=0;}


    // Returns the chisq at current basis and B vector:
    double getChisq() const {return fitIsValid ? chisq : -1.;}
    int    getDOF() const {return fitIsValid ? DOF : -1; }

    const  LVector& getB() const {return b;}
    // Get the (SVD) of the inverse covariance matrix for the b coefficients.
    // F = V^T S V, V is unitary.
    tmv::ConstVectorView<double> fisherSV_S() const;
    tmv::ConstMatrixView<double> fisherSV_V() const;
    // Get covariance matrix for b (inverse of F):
    tmv::SymMatrix<double> covar();

    // Report 1st coeff and its error if it were only fitted term:
    void b00(double& value, double& var) const;
    // Get flux significance, using monopole terms up to maxP - at fixed basis
    double significance(int maxP=-1);
    // Covariance of E's, marginalizing over other varying basis elements
    void   shapeVariance(double& sigE1, double& sigE2, double& covar) ;
    // Error on size, marginalizing over other basis elements
    double muError();

    // Roundness test execution - returns elements currently set to vary
    DVector Mtests();
    DMatrix dMdE(double xyscale=1.);  // optional scale factor to apply to x0/y0 shifts

  private:
    // Parameters defining what is desired from this fit:
    int order;
    Ellipse basis;
    Ellipse maskBasis;

    bool centering;
    bool shearing;
    bool dilating;

    // The data, and the function class used to fit data:
    typename std::list<FitExposure<T> > exposures;

    // Data arrays:
    DVector* x;
    DVector* y;
    DVector* dsig;
    DVector* invsig;

    // Fitting results:
    int     flags;
    LVector b;
    double  chisq;
    int     DOF;
    double  f00;	// First coeff if it were the only thing being fit
    double  var00;      // ...and its variance.
    // The design matrix:
    DMatrix* A;

    // Roundness tests and their derivative matrices:
    DMatrix* M;
    vector< DMatrix*> MG;
    void updateMG();
    void updateM();
    // Return a vector listing currently active roundness tests
    vector<LVector::GType> buildTestIndices();

    /////////////////////
    // State indicators & manipulation
    void setDefaults();

    // Are masks chosen and data acquired
    bool    dataAreValid;
    void    invalidateData();
    bool    acquireData();

    // Are b vector, alpha, chisq, etc., valid for current data & basis & order?
    bool    fitIsValid;
    void    invalidateFit();
    
    // Does current basis b, alpha, etc. represent solution for round/shear/size ?
    bool    isSolved;

    // Magic numbers for fitting process:
    double tolerance;
    double maxMaskMu;
    int    maxIterations;

    // How large is fitting region (in units of maskBasis)
    double maskSigma;

    // See how large/small currently optimal mask is compared
    // to the current actual mask. 
    void maskRatios(double& smallest, double& largest);

    // And a we'll keep a vector of them for each exposure
    vector<CrudeMap> maps;
  };		 
  
} // namespace
#endif // GLSIMPLE_H
