// $Id: FDNT.h,v 1.24 2012/01/03 15:43:52 dgru Exp $

// Fourier-Domain Null Test shape measurement codes

#ifndef FDNT_H
#define FDNT_H

#include <cmath>
#include <list>

#include "Std.h"
#include "UseTMV.h"
#include "TMV_SymBand.h"
#include "Shear.h"
#include "Image.h"
#include "Laguerre.h"
#include "SBProfile.h"
#include "GTable.h"
#include "GLSimple.h"
#include "Lance.h"
#include "Random.h"
#include "FitFlags.h"

namespace laguerre {

  class FDNTError: public std::runtime_error {
  public:
    FDNTError(const string &m=""): std::runtime_error("FDNT Error: " +m) {}
  };

  // Forward declaration
  template <class T>
  class FDNT;

  // Interface class to circular functions
  // Pass arguments as square of k to avoid sqrts.
  class CircleK {
  public:
    CircleK(): scaleFactorSq(1.) {}
    virtual double operator()(double ksq) const=0;
    virtual double dFdlnScale(double ksq) const {
      const double DLNS=0.01;
      return ((*this)(ksq*scaleFactorSq*(1+DLNS))
	      -(*this)(ksq*scaleFactorSq*(1-DLNS)) ) / DLNS;
    }
    void setScaleFactor(double s) {scaleFactorSq=s*s;}
    double getScaleFactor() const {return sqrt(scaleFactorSq);}
    virtual ~CircleK() {}
  protected:
    double scaleFactorSq;
  };

  // A Gaussian version of above:
  class GaussKsq: public CircleK {
  public:
    GaussKsq(double sigma=1.) {setScaleFactor(sigma);} // DG: are you sure the argument of setScaleFactor should be sigma and not sigma^2
    double operator()(double ksq) const {return exp(-0.5*ksq*scaleFactorSq);}
    double dFdlnScale(double ksq) const {
      return -ksq*scaleFactorSq*exp(-0.5*ksq*scaleFactorSq);
    }
  };

  // ...and one that is Fourier Transform of exponential disk:
  class ExpDiskKsq: public CircleK {
  public:
    ExpDiskKsq(double r0=1.) {setScaleFactor(r0*r0);}
    double operator()(double ksq) const {return pow(1+ksq*scaleFactorSq, -2.5);}
    double dFdlnScale(double ksq) const {
      return -5.*ksq*scaleFactorSq*pow(1+ksq*scaleFactorSq, -3.5);
    }
  };

  // A circular function that azimuthally averages an SBProfile
  // and SQUARES it (for the T_0^2 optimal weight function):
  // Will build a lookup table for spline interpolation
  // and also handles nulls in MTF to some extent.
  class CircularizedSB: public CircleK {
  public:
    CircularizedSB(const sbp::SBProfile& sb); 
    double operator()(double ksq) const {
      if (ksq*scaleFactorSq >= maxksq) return 0.;
      return tab(ksq*scaleFactorSq);
    }
  private:
    double maxksq;
    GTable<> tab;
  };

  // Structure to give all the info we want about a PSF:
  class PSFInformation {
  public:
    PSFInformation(const sbp::SBProfile& T_, Ellipse e_, CircleK* T0sq_):
      T(T_), shape(e_), T0sq(T0sq_), ownT0sq(false) {}
    PSFInformation(const sbp::SBProfile& T_, Ellipse e_):
      T(T_), shape(e_), T0sq(new CircularizedSB(T_)), ownT0sq(true) {}
    ~PSFInformation() {if (ownT0sq) delete T0sq;}

    // The PSF itself (transfer function):
    const sbp::SBProfile& T;
    // Rough description of size and shape of PSF:
    Ellipse shape;
    // Circularized function to use for weight, nominally inverse-square of the PSF:
    CircleK* T0sq;
  private:
    bool ownT0sq;
    // Hide copy constructor:
    PSFInformation(PSFInformation& rhs): T(rhs.T) {}
  };

  struct TestArrays;	// scratch vectors, defined in FDNT.cpp

  template <class T=float>
  class ExposureGroup {
    friend class FDNT<T>;
  protected:
    ExposureGroup(PSFInformation& psf_);
    list<FitExposure<T> > felist;
    PSFInformation& psf;
    Position<double> phaseCenter;	// world coords at origin of Fourier transform
  public:
    // Standard definitions of indices for ellipse elements
    // within the test arrays (iE1,iE2,iSize,iX,iY).
    // nRe and nIm are number of real (iE1,iE2,iSize) and imag (iX,iY) tests in k space.
    // nTests is total number of tests.
    static void setIndices(bool centerIsFixed,
			   bool sizeIsFixed,
			   int& iE1, int& iE2,
			   int& iSize,
			   int& iX, int& iY,
			   int& nRe, int& nIm,
			   int& nTests);
    // Return null tests and calculate quantities of interest.
    // Assumes properly set scale factors in wg and T0sq going in.
    DVector tests(Shear targetS,		// Shear to test
		  const CircleK& wg,		// Galaxy weight function
		  DMatrix* dMdE,		// Deriv w.r.t. ellipse params, if !=0
		  DVector* dMdlns,		// Derivs w.r.t. T0 shrink factor, if !=0
		  tmv::SymMatrix<double>* covM,	// Covariance matrix, if !=0
		  bool centerIsFixed,		// true to hold centroid fixed
		  bool sizeIsFixed,		// true to hold wg weight scale fixed.
		  bool exactdMdE);		// false to limit derivs to m<=2 multipoles
    virtual ~ExposureGroup();
    // Do all preparatory calculations and return
    // weight to apply to null tests of this ExposureGroup:
    virtual double prepare(Ellipse startBasis) =0;

    // Rotate phases to define new world-coord phase center
    virtual void setPhaseCenter(double x0, double y0)=0;
    Position<double> getPhaseCenter() const {return phaseCenter;}

    // Get FDNT-weighted flux and its variance
    void wtFlux(Shear targetS,			//Shear to test
		const CircleK& wg,		//Galaxy weight function
		double& f, double& varf);

    int getFlags() const {return flags;}

    // For investigation purposes: draw the real-space version of the e1 test
    // worldBasis used to choose pixels to draw on.
    Image<> drawFilter(Shear targetS,
		       const CircleK& wg,
		       Ellipse worldBasis,
		       bool ise1);

  protected:

    // Flags from processing this exposure:
    int flags;

    // k vectors
    DVector* ktruex;
    DVector* ktruey;

    TestArrays* arrayPtr;

    // Real & imag parts of FT of deconvolved object
    DVector* deconvRe;
    DVector* deconvIm;
    // Factor applied to each FT element to make an integral over d^2k:
    double ftnorm;
    // Index of DC element, which is has no conjugate term
    // and hence doesn't get 2x in integals:
    // [note that other self-conjugate k's are at edges of k domain, where they
    // should not be getting significant weight.
    int dcIndex;

    // Derivs of FT w.r.t. translation:
    DMatrix* dRedxy;
    DMatrix* dImdxy;
    void createdXdY();

    // Calculate the covariance matrix of union of two vectors, one is
    // mRe * deconvRe  - where mRe dimensions are (nRe x nk), then
    // mIm * deconvIm  - (nIm x nk).
    // Return a symmetric cov matrix of dimension nRe + nIm.
    virtual tmv::SymMatrix<double> covarianceOf(const DMatrix& mRe, 
						const DMatrix& mIm)=0;

  };

  // Exposure group with transform to be executed via GL decomposition
  template <class T=float>
  class ExposureGroupGL: public ExposureGroup<T> {
  private:
    LVector b;	// tells us the order
    double maskSigma;

    // (square root of) covariance matrices for Re/Im Fourier
    // coefficients
    DMatrix* bReV;
    DMatrix* bImV;

    using ExposureGroup<T>::ktruex;
    using ExposureGroup<T>::ktruey;
    using ExposureGroup<T>::deconvRe;
    using ExposureGroup<T>::deconvIm;
    using ExposureGroup<T>::ftnorm;
    using ExposureGroup<T>::dcIndex;

  public:
    ExposureGroupGL(PSFInformation& psf_,
		    int order,
		    double maskSigma_=-1.);
    ~ExposureGroupGL();
    virtual double prepare(Ellipse startBasis);

    void setMaskSigma(double maskSigma_) {maskSigma=maskSigma_;}
    virtual tmv::SymMatrix<double> covarianceOf(const DMatrix& mRe, 
						const DMatrix& mIm);
    virtual void setPhaseCenter(double x0, double y0);
  };

  // Exposure group to use Fourier Transform.  
  // Not using multiple exposures yet.
  template <class T=float>
  class ExposureGroupFT: public ExposureGroup<T> {
  private:
    // Variance of real transform elements
    tmv::DiagMatrix<double>* varRe;

    using ExposureGroup<T>::ktruex;
    using ExposureGroup<T>::ktruey;
    using ExposureGroup<T>::deconvRe;
    using ExposureGroup<T>::deconvIm;
    using ExposureGroup<T>::ftnorm;
    using ExposureGroup<T>::dcIndex;

  public:
    ExposureGroupFT(PSFInformation& psf_);
    ~ExposureGroupFT();
    virtual double prepare(Ellipse startBasis);
    virtual tmv::SymMatrix<double> covarianceOf(const DMatrix& mRe, 
						const DMatrix& mIm);
    virtual void setPhaseCenter(double x0, double y0);
  };

  // Now the class that does the FDNT tests
  template <class T=float>
  class FDNT {
  public:
    // Constructors: nativeBasis is estimated mean image size/posn
    // order gives GL order, or <=0 means use direct FT
    FDNT(list<FitExposure<T> > felist,
	 vector<PSFInformation*> psfv,
	 Ellipse nativeBasis_,
	 int order=12);
    // Constructor for one exposure, dressed up
    FDNT(FitExposure<T> fe,
	 PSFInformation& psf,
	 Ellipse nativeBasis_,
	 int order=12);
    // Simplest constructor: trivial coordinate map, and we are
    // responsible for making T0 & measuring PSF
    FDNT(const Image<T> sci, const Image<T> invvar, 
	 const sbp::SBProfile& psf, 
	 Ellipse nativeBasis_,
	 int order=12);
    ~FDNT();

    // Update nativeBasis from current value by doing a joint
    // GL decomposition of all native images - want this to 
    // get a common centroid.
    Ellipse GLAll(int order=4);

    // Acquire all data and decompositions from all images
    bool prepare();

    // Hold centroid fixed during searches:
    void setCentering(bool rectr) {centerIsFixed=!rectr;}
    bool getCentering() const {return !centerIsFixed;}
    void setSizing(bool resize) {sizeIsFixed=!resize;}
    bool getSizing() const {return !sizeIsFixed;}

    // Calculate probability of null being within dEta1 dEta2 of this shear
    // and the cov matrix for shear at this point.
    // Sets CentroidMismatch and/or SizeMismatch flags if max likelihood
    //  positions are far from nominal point (linear approx may fail...)
    // Marginalization over centroid/size decided by current flags.
    double logProbability(Shear targetS,	
			  tmv::SymMatrix<double>& covE);
    // In this one, the marginalization over ln(sigma) is done numerically
    // if needed.
    double logProbability2(Shear targetS,	
			   double& fractionTooSmall,
			   double& fractionTooBig,
			   tmv::SymMatrix<double>& covE);

    // Maximum likelihood search, starting at current basis:
    // check flags for success, warnings, failures
    Shear shape(double& logLikelihood,
		tmv::SymMatrix<double>& covE);

    // Use Lance Miller's algorithm to map likelihood:
    Shear shape2(double& logLikelihood,
		tmv::SymMatrix<double>& covE);

    /***** Worry about these later ??? 
    // Return value of E1/E2 null tests (simple sum over ExposureGroups)
    // If centering and resizing are on, returns values at size/center
    // that minimize norm of E1/E2.  Otherwise gives E1/E2
    // Optional derivative matrix.
    DVector testE1E2(Shear targetS, DMatrix* dMdE);

    // Dogleg search for a null m1/m2, starting at current basis
    // Use simple sum of tests instead of trying to weight.
    bool doglegNull(double etaTol=0.05);
    *****/

    // Get flux under the FDNT filter, and its variance
    // Exposure groups weighted by inverse flux variance.
    void wtFlux(double& f,
		double& varf);

    // Get responsivity bias due to ellipticity gradients
    double shrinkResponse(Shear targetS);

    int getFlags() const {return flags;}
    void setFlag(int f) { flags = flags | f; return;}
    // Return current basis in use for galaxy weighting
    Ellipse getBasis() const {return galaxyBasis;}
    // Set the size used for galaxy weighting:
    void setIntrinsicSize(double sigma) {galaxyBasis.setMu(log(sigma));}
    // Fix an estimate for galaxy center (which will force
    // re-phasing of FT's to be centered at this position)
    void setCenter(double x0, double y0);
    // Find out where origin of FT's is currently set:
    Position<double> getPhaseCenter() const;

    void setMaskSigma(double sig);

    static double shrinkFactor(Shear targetS) ;
    static double dlnShrinkdEta(Shear targetS);

    // Draw the ellipticity filter in real space.
    // worldBasis describes general size/location of the
    // object for purpose of deciding what pixel range of postage stamp
    // should be.
    Image<> drawFilter1(Shear targetS, Ellipse worldBasis) {
      //cout << "drawing filter in " << this << " exposure group " << veg.front() << endl;
      return veg.front()->drawFilter(targetS, *wg, worldBasis, true);
    }
    Image<> drawFilter2(Shear targetS, Ellipse worldBasis) {
      return veg.front()->drawFilter(targetS, *wg, worldBasis, false);
    }

    // Some things to peek into algorithm performance:
    int getEvaluationCount() const {return evaluationCount;}
    void resetEvaluationCount() {evaluationCount=0;}
    int getETrialCount() const {return eTrialCount;}
    double getSampleDensity() const {return sampleDensity;}
    double getTotalProbability() const {return totalProbability;}

  private:
    static ran::UniformDeviate u;
    int flags;
    int evaluationCount;
    int eTrialCount;
    double sampleDensity;
    double totalProbability;

    vector<ExposureGroup<T>*> veg;
    vector<bool> useEG;   // Set to true if ExposureGroup has valid tests.

    CircleK* wg;	// galaxy weight function
    bool centerIsFixed;
    bool sizeIsFixed;
    // Pointers to classes we might own & need to destroy:
    PSFInformation* localPSFInfo;
    FitExposure<T>* localFitExposure;

    double maskSigma;	// Default mask size for GL fits
    double leastPSFSigma;

    Ellipse nativeBasis; // rough observed world location / size of galaxy
    Ellipse galaxyBasis; // estimated intrinsic shape & center of galaxy, and weight size

    void setShrinkFactor(double s);

    // Internal routine to sum over ExposureGroups
    // Produces Fisher as total Fisher matrix over E variables;
    //   DCinvT is sum of (dt/dE)^T C_t^{-1} t  where t is test vector
    //   DCinvdTdlnS is sum (dt/dE)^T C_t^{-1} (dt/d(ln s)).
    // Passing in zero pointers bypasses evaluation.
    // Size of weight function taken from galaxyBasis.  phaseCenter
    //  is origin of tests, unless centerIsFixed, in which case
    //  the phaseCenter is reset to the galaxyBasis center.
    void  sumTests(Shear targetS,
		   tmv::SymMatrix<double>& Fisher,
		   DVector* DCinvT,
		   DVector* DCinvdTdlnS);

    // Sum tests and return test values & covariance matrix marginalized over
    //  centroid position, if !centerIsFixed.  Outputs are either nRe= 2 or 3 dimensional
    //  depending on whether sizeIsFixed.
    // If marginalizing, centroid of galaxyBias is updated to new max likelihood posn
    void marginalizedCentroid(Shear targetS,
			      DVector& dE23,
			      tmv::SymMatrix<double>& covE23);

    // Sum tests, marginalize over both centroid and/or size as specified by current flags.
    // Outputs are test values and covariance matrix over E1/E2.
    // Updates galaxyBasis to have size and center of max likelihood point, if marginalizing
    void marginalizedSize(Shear targetS,
			  DVector& dE2,
			  tmv::SymMatrix<double>& covE2);

    void processSample(Sample& s);

    void setDefaults();

  };
    
} // namespace laguerre
#endif   // FDNT_H
