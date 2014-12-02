// $Id: FDNT.cpp,v 1.21 2012/01/03 15:40:09 dgru Exp $

#include "FDNT.h"
#include "FitFlags.h"
#include "Brent.h"
#include "Solve.h"
#include "Random.h"
#include <set>

// *** #define DEBUG

using namespace laguerre;

// Flag if peak-likelihood centroid is more than this many native sigmas
// away from phase center (i.e. linear dependence on center may fail)
const double CENTROID_MISMATCH_THRESHOLD=0.3;
// Flag if ln(galaxy sigma) which maximizes likelihood is more than this far
// from the size used to construct tests:


const double SIZE_MISMATCH_THRESHOLD=0.1;
// Largest allowed motion in ln(galaxy sigma)
const double MAXIMUM_SIZE_SHIFT=2*SIZE_MISMATCH_THRESHOLD;
const double MAXIMUM_SIZE_FACTOR = 2.;
const double MINIMUM_SIZE_FACTOR = 0.2;

const int MAX_ITERATIONS = 40;
const double MINIMUM_TRUST_RADIUS=0.03;
const double INITIAL_TRUST_RADIUS=0.3;
// Will not search for ellipticity solution beyond this eta:
//const double MAXIMUM_ETA=2.8;
const double MAXIMUM_DELTA_ETA=1.5;
// More than this probability beyond this eta sets flag:
//const double WARN_ETA_THRESHOLD=2.5;
//const double WARN_ETA_PROBABILITY=0.5;

// (log of) factor below PSF size that object is considered unresolved
const double TOO_SMALL = log(2.);
// (log of) factor above native size considered too large
const double TOO_BIG = log(1.5);

// Fraction of likelihood with flag set that 
// will trigger output flag
const double CENTROID_FLAG_THRESHOLD = 0.2;
// Or fraction with galaxy sigma outside sensible bounds
const double SIZE_PILEUP_THRESHOLD = 0.2;

template <class T>
void
FDNT<T>::setDefaults() {
  wg = 0;
  localPSFInfo = 0;
  localFitExposure = 0;
  centerIsFixed = false;
  sizeIsFixed = true;
  maskSigma = -1.;
  flags = 0.;
  evaluationCount = 0;
}

template <class T>
FDNT<T>::~FDNT() {
  if (wg) delete wg;
  if (localPSFInfo) delete localPSFInfo;
  if (localFitExposure) delete localFitExposure;
  for (int i=0; i<veg.size(); i++) delete veg[i];
}

template <class T>
double
FDNT<T>::shrinkFactor(Shear targetS) {
  return exp(0.5*targetS.getEta());
}

template <class T>
double
FDNT<T>::dlnShrinkdEta(Shear targetS) {
  return 0.5;
}

template <class T>
void
FDNT<T>::setMaskSigma(double sig) {
  maskSigma = sig;
  for (int i=0; i<veg.size(); i++) {
    ExposureGroupGL<T>* ptr = dynamic_cast<ExposureGroupGL<T>*> (veg[i]);
    if (ptr) ptr->setMaskSigma(sig);
  }
}

template <class T>
void
FDNT<T>::setCenter(double x0, double y0) {
  galaxyBasis.setX0(Position<double>(x0,y0));
  // Reset the native basis too, in case we have not yet run prepare()
  nativeBasis.setX0(Position<double>(x0,y0));
  for (int i=0; i<veg.size(); i++)
    veg[i]->setPhaseCenter(x0,y0);
}

template <class T>
Position<double>
FDNT<T>::getPhaseCenter() const {
  return veg.front()->getPhaseCenter();
}


// Constructor with single image
template <class T>
FDNT<T>::FDNT(const Image<T> sci, const Image<T> invvar, 
	      const sbp::SBProfile& psf, 
	   Ellipse nativeBasis_,
	   int order): nativeBasis(nativeBasis_) {

  setDefaults();
  localFitExposure = new FitExposure<T>(sci, invvar, 0);

  // Measure the PSF using a 4th-order GL fit in 4-sigma mask.
  Image<> psfImg=psf.draw();
  double dx;
  Image<> wt = psfImg.duplicate();
  wt = 1.;
  Ellipse psfBasis(0., 0., log( (psfImg.XMax()-psfImg.XMin())/10.), 0., 0.);
  GLSimple<T> gl(psfImg, wt, psfBasis, 4);
  gl.setMaskSigma(4.);
  gl.solve();
  // Look for flags other than BDegeneracy, which will get set
  // since we're not very careful about noise level
  if (gl.getFlags() & ~BDegeneracy) 
    FormatAndThrow<FDNTError>() << "GL flags set while measuring PSF: " 
				<< gl.getFlags();
  psfBasis = gl.getBasis();
  if (!psfImg.getHdrValue("DX",dx))
    throw FDNTError("Missing DX while measuring PSF size");
  psfBasis.setMu( psfBasis.getMu() + log(dx));
  localPSFInfo = new PSFInformation(psf, psfBasis);

  if (order >0) {
    veg.push_back(new ExposureGroupGL<T>(*localPSFInfo, order, maskSigma));
  } else {
    veg.push_back(new ExposureGroupFT<T>(*localPSFInfo));
  }
  veg.back()->felist.push_back(*localFitExposure);
  useEG.push_back(true);
}

template <class T>
FDNT<T>::FDNT(FitExposure<T> fe,
	      PSFInformation& psf,
	      Ellipse nativeBasis_,
	      int order): nativeBasis(nativeBasis_) {
  setDefaults();
  if (order > 0)
    veg.push_back( new ExposureGroupGL<T>(psf, order, maskSigma) );
  else
    veg.push_back( new ExposureGroupFT<T>(psf));
  
  veg.back()->felist.push_back(fe);
  useEG.push_back(true);
}

template <class T>
FDNT<T>::FDNT(list<FitExposure<T> > felist,
	      vector<PSFInformation*> psfv,
	      Ellipse nativeBasis_,
	      int order): nativeBasis(nativeBasis_) {
  setDefaults();
  
  // Put each FitExposure into an ExposureGroup that uses
  // its PSF.
  for ( typename list<FitExposure<T> >::iterator ife = felist.begin(); 	// for all exposures do
	ife != felist.end();
	++ife) {
    int currentPSFCode=ife->psfCode;					// create new ExposureGroup if this exposure has new PSF
    if (currentPSFCode >= psfv.size())
      FormatAndThrow<FDNTError>() << "Requested psfCode " << currentPSFCode
				  << " is not in PSFInfo vector given to constructor";
    typename vector<ExposureGroup<T>* >::iterator ieg;
    for (ieg=veg.begin(); ieg != veg.end(); ++ieg)
      if ((*ieg)->felist.front().psfCode == currentPSFCode) break;
    
    if (ieg==veg.end()) {
      // New exposure group: 
      if (order > 0)
	veg.push_back( new ExposureGroupGL<T>(*psfv[currentPSFCode], order, maskSigma) );
      else
	veg.push_back( new ExposureGroupFT<T>(*psfv[currentPSFCode]));

      veg.back()->felist.push_back(*ife);
      useEG.push_back(true);
    } else {
      // Add this FitExposure to existing ExposureGroup:
      (*ieg)->felist.push_back(*ife);
    }
  }
  Assert(useEG.size()==veg.size());
}

template <class T>
Ellipse
FDNT<T>::GLAll(int order) {
  // Collect all our exposures into a common fitting list:
  list<FitExposure<T> > allfe;
  for (int i=0; i<veg.size(); i++) {
    for (typename list<FitExposure<T> >::iterator ife=veg[i]->felist.begin();
	 ife != veg[i]->felist.end();
	 ++ife)
      allfe.push_back(*ife);
  }
  // Execute GL fitting at requested order
  GLSimple<T> gl(allfe, nativeBasis, order);
  gl.setMaskSigma(maskSigma);
  // Fix the origin if it's not free to move:
  if (centerIsFixed) gl.setCentering(false);

  if (!gl.solve()) {
    // Fatal error in fitting
    flags = GLFailure + gl.getFlags();
    return Ellipse();
  } else {
    // Some valid shape produced:
    flags = gl.getFlags();
    nativeBasis = gl.getBasis();
    return nativeBasis;
  }
}

// Prepare all images.  Returns true if there is successful preparation
// of *any* image.
template <class T>
bool
FDNT<T>::prepare() {
  // Get rid of any previous preparation results
  flags = 0;  // clear any previous problems
  bool success = false;
  double sumw = 0.;
  double sumpsfxx=0., sumpsfyy=0., sumpsfxy=0.;
  leastPSFSigma = 1.e20;
  for (int i=0; i<veg.size(); i++) {
    double sn = veg[i]->prepare(nativeBasis);
    // Edge flags propagate into our flags:
    flags += veg[i]->getFlags() & (Edge);
    if (sn==0.) {
      useEG[i] = false;  // Do not use this ExposureGroup any further.
    } else {
      useEG[i] = true;
      success = true;
      double wt = sn*sn;	// approx PSF estimate weighted by S/N ratios
      sumw += wt;
      leastPSFSigma = MIN(leastPSFSigma, exp(veg[i]->psf.shape.getMu()));
      double xxplusyy = 2.*exp(2.*veg[i]->psf.shape.getMu())
	* cosh(veg[i]->psf.shape.getS().getEta());
      sumpsfxx += wt * xxplusyy * (1+veg[i]->psf.shape.getS().getE1());
      sumpsfyy += wt * xxplusyy * (1-veg[i]->psf.shape.getS().getE1());
      sumpsfxy += wt * xxplusyy * veg[i]->psf.shape.getS().getE2() /2.;
    }
  }
  // If no exposures succeeded, report the union of all flags:
  if (!success)
    for (int i=0; i<veg.size(); i++) 
      flags |= veg[i]->getFlags();

  if (success) {
    // Choose an initial size for galaxy weight function
    // ??? This could clearly use some refinement
    // ??? Need a shear-independent assignment; subtract moments?
    if (wg) delete wg;
    double galaxySigma; // = exp(2.*nativeBasis.getMu()) - psfsq;
    {
      // Calculate galaxy size via subtraction of squares
      double xxplusyy = 2.*exp(2.*nativeBasis.getMu())
	* cosh(nativeBasis.getS().getEta());
      double ixx = xxplusyy * (1+nativeBasis.getS().getE1()) - sumpsfxx/sumw;
      double iyy = xxplusyy * (1-nativeBasis.getS().getE1()) - sumpsfyy/sumw;
      double ixy = xxplusyy * nativeBasis.getS().getE2()/2. - sumpsfxy/sumw;
      galaxySigma = sqrt( MAX(0., ixx*iyy - 4*ixy*ixy))/2.;
    }
    galaxySigma = sqrt(galaxySigma);
    galaxySigma = MAX( galaxySigma, MINIMUM_SIZE_FACTOR*leastPSFSigma);
    // ??? turn off sigma iteration for small sigmas ???
    //**/cerr << "galaxy sigma " << galaxySigma << endl;
    wg = new GaussKsq(galaxySigma);
    //**/ ???wg = new ExpDiskKsq(galaxySigma);
    galaxyBasis = nativeBasis;
    galaxyBasis.setMu(log(galaxySigma));
  }
  return success;
}

template <class T>
void
FDNT<T>::setShrinkFactor(double s) {
  for (int i=0; i<veg.size(); i++)
    veg[i]->psf.T0sq->setScaleFactor(s);
}

template <class T>
void
FDNT<T>::sumTests(Shear targetS,
		  tmv::SymMatrix<double>& Fisher,
		  DVector* DCinvT,
		  DVector* DCinvdTdlnS) {
  
  //cout << "sumTests called with " << targetS << Fisher << DCinvT << DCinvdTdlnS << endl;
  
  // Get dimensions
  int iE1, iE2, iSize, iX, iY, nRe, nIm, nTests;
  ExposureGroup<float>::setIndices(centerIsFixed, sizeIsFixed,
				   iE1, iE2, iSize, iX, iY, nRe, nIm, nTests);

  // Check input vector sizes
  Assert(Fisher.nrows()==nTests);
  Fisher.setZero();
  if (DCinvT) {
    Assert(DCinvT->size()==nTests);
    DCinvT->setZero();
  }
  if (DCinvdTdlnS) {
    Assert(DCinvdTdlnS->size()==nTests);
    DCinvdTdlnS->setZero();
  }
  
  Assert(wg);
  // Set wg scale from galaxyBasis
  wg->setScaleFactor( exp(galaxyBasis.getMu()) );
  //**/cerr << "*** logProb using sigma= " << wg->getScaleFactor() << endl;
  // Choose shrink Factor, apply to all T0sq filters:
  setShrinkFactor( shrinkFactor(targetS) );

  // If we are not marginalizing over centroid, need to
  // fix center of tests to current position
  if (centerIsFixed) setCenter(galaxyBasis.getX0().x,
			       galaxyBasis.getX0().y);


  // Vectors/matrices for individual ExposureGroups
  tmv::SymMatrix<double> Ceg(nTests, 0.);
  DMatrix dtdE(nTests, nTests, 0.);
  DVector dTdlnS(nTests, 0.);
  DVector t(nTests, 0.);
  
  evaluationCount++; 

  // Try block for catching all matrix singularities:
  try {
    for (int i=0; i<veg.size(); i++) {
      //cout << "// For each exposure group:" << endl;
      if (!useEG[i]) continue;

      //cout << "//   Get raw test data" << endl;
      t = veg[i]->tests(targetS, *wg,
			&dtdE,
			DCinvdTdlnS ? &dTdlnS : 0 ,
			&Ceg,
			centerIsFixed,
			sizeIsFixed,
			true);	// ??? exact dMdE in use always
      //cout << "//   Make D/C, make Fisher and dE summands" <<  endl;
      DMatrix CinvD = dtdE / Ceg;
      
      //cout << "//   accumulate into full sum set" << endl;
      //   ...following line does Ftot += dtdE.transpose * CinvD
      tmv::SymMultMM<true>(1., dtdE.transpose(), CinvD.view(), Fisher.view());
      if (DCinvT) (*DCinvT) += t * CinvD;
      if (DCinvdTdlnS) (*DCinvdTdlnS) += dTdlnS * CinvD;

      if (!(Fisher(0,0)==Fisher(0,0))) {
	cerr << "NaN Fisher.  Last Ceg:\n" << Ceg
	     << " dtdE\n " << dtdE
	     << " t\n" << t
	     << " galaxyBasis " << galaxyBasis
	     << endl;
	flags |= Singularity;
	throw FDNTError("NaN Fisher in sumTests()"); 
      }
    }
  } catch (tmv::Singular) {
    flags |= Singularity;
    throw FDNTError("Singular matrix in sumTests()");
  }
}

template <class T>
void
FDNT<T>::marginalizedCentroid(Shear targetS,
			      DVector& dE23,
			      tmv::SymMatrix<double>& covE23) {
  int iE1, iE2, iSize, iX, iY, nRe, nIm, nTests;
  ExposureGroup<float>::setIndices(centerIsFixed, sizeIsFixed,
				   iE1, iE2, iSize, iX, iY, nRe, nIm, nTests);
				   
  //cout << "marginalizedCentroid called: " << targetS << " " << dE23 << " " << covE23 << endl;
  // Check input vector sizes
  Assert(dE23.size()==nRe);
  Assert(covE23.nrows()==nRe);

  // Get the full-dimensional information
  tmv::SymMatrix<double> Fisher(nTests,0.);
  DVector DCinvT(nTests,0.);
 
  
  sumTests(targetS,
	   Fisher,
	   &DCinvT,
	   0);
  try {
    Fisher.saveDiv();
    tmv::SymMatrix<double> Ctot = Fisher.inverse(); 
    DVector dE = DCinvT / Fisher;
    //  If we are marginalizing over centroid, set flag if
    // we end up too far from starting phase center
    if (iX>=0) {
      // ??? Is this ok ???
      DVector xyShift = dE.subVector(iX, iY+1)
	+ (Fisher.subMatrix(iX,iY+1,0,nRe) / Fisher.subSymMatrix(iX, iY+1) ) * dE.subVector(0, nRe);
      if (hypot(xyShift[0],xyShift[1]) > CENTROID_MISMATCH_THRESHOLD * exp(nativeBasis.getMu())) {
	/*/
	cerr << "CentroidMismatch set #1.  Difference is: "
	     << hypot(xyShift[0],xyShift[1]) - CENTROID_MISMATCH_THRESHOLD * exp(nativeBasis.getMu())
	     << endl;
	/*/
	flags |= CentroidMismatch;
      } else {
	flags &= ~CentroidMismatch;
      }
      // Put likelihood-maximizing centroid into galaxyBasis:
      galaxyBasis.setX0( getPhaseCenter() - Position<double>(xyShift[0],xyShift[1]) );
    }

    //  Marginalize over centroid & size by extracting just the E1/E2/Size
    // parts of covariance and tests.
    covE23 = Ctot.subSymMatrix(0, nRe);
    dE23 = dE.subVector(0, nRe);
  } catch (tmv::Singular) {
    flags |= Singularity;
    throw FDNTError("Singular Fisher matrix in marginalizedCentroid()");
  }

}


template <class T>
double
FDNT<T>::shrinkResponse(Shear targetS) {
  if (targetS.getE()==0.) return 0.;	// no bias at circular test

  // Get dimensions, set up full Fisher and dE matrices
  int iE1, iE2, iSize, iX, iY, nRe, nIm, nTests;
  ExposureGroup<float>::setIndices(centerIsFixed, sizeIsFixed,
				   iE1, iE2, iSize, iX, iY, nRe, nIm, nTests);

  tmv::SymMatrix<double> Ftot(nTests, 0.);
  DVector DCinvdTdlnS(nTests, 0.);
  sumTests(targetS, Ftot, 0, &DCinvdTdlnS);
  try {
    DVector dEdlns = DCinvdTdlnS / Ftot;
    // Take component along radius in e plane
    double d = dEdlns[iE1] * targetS.getE1() + dEdlns[iE2]*targetS.getE2();
    return d * dlnShrinkdEta(targetS) / targetS.getE();
  } catch (tmv::Singular) {
    flags |= Singularity;
    throw FDNTError("Singular matrix in shrinkResponse()");
  }
}

template <class T>
void
FDNT<T>::processSample(Sample& s) {
  tmv::SymMatrix<double> covE(2);
  double eta1, eta2;
  s.getXY(eta1,eta2);
  //cout << "processSample: gotXY" << endl;
  Shear trialS;
  trialS.setEta1Eta2(eta1,eta2);
  //cout << "processSample: calling logProbability2" << endl;
  s.lnProb = logProbability2(trialS, s.fractionTooSmall, s.fractionTooBig, covE);
  //cout << "processSample: gettingMu" << endl;
  s.sigma = galaxyBasis.getMu();
  s.centroidFlag = CentroidMismatch & flags; 
    // logProbability2 called marginalizedCentroid, which re-sets CentroidMismatch flag to flag for this sample
    // need to save it in the sample though, since the global flag will be overwritten by the next sample
  s.covE = covE;
}

template <class T>
ran::UniformDeviate
FDNT<T>::u;

template <class T>
Shear
FDNT<T>::shape2(double& logLikelihood,
		tmv::SymMatrix<double>& covE) {
  //cout << "shape2 starting up" << endl;
  //const double MINSTEP=0.01;
  //const int START_LEVEL=5;	// initial grid 2^5 * MINSTEP = 0.32
  const double MINSTEP=8e-5;
  const int START_LEVEL=12;	// initial grid 2^12 * MINSTEP = 0.32768
  //**/const int MINPOINTS=30;
  const int MINPOINTS=50;
  //**/  const double PROB_THRESHOLD=-log(0.01);
  const double PROB_THRESHOLD=-log(0.002);

  const double lnSigmaFloor = log(leastPSFSigma)-TOO_SMALL;
  const double lnSigmaCeiling = nativeBasis.getMu() + TOO_BIG;

  // Create origin node - precisely at shape of current galaxyBasis
  double eta1o, eta2o;
  galaxyBasis.getS().getEta1Eta2(eta1o,eta2o);
  Sample origin(eta1o, eta2o, MINSTEP, START_LEVEL);
  processSample(origin); 
  //cout << "shape2 did processSample" << endl;

  // Set max probability (actually min of negative log prob)
  double minProb = origin.lnProb;

  // Create current list
  list<Sample> current;
  set<Sample> tried;
  current.push_back(origin);
  tried.insert(origin);

  int radius = 0;
  bool gotOne;
  do {
    // Look at all neighbors of current valid points
    // Keep growing until no new points that are in bounds and
    // above threshold.
    gotOne = false;
    for (list<Sample>::iterator i=current.begin(); i!=current.end(); ++i) {
      // Skip if below threshold
      if (i->lnProb >= minProb + PROB_THRESHOLD) continue;
      // Get all neighbors of this point
      list<Sample> n=i->neighbors();
      for (list<Sample>::iterator j=n.begin(); j!=n.end(); ++j) {
	double x,y;
	j->getXY(x,y);
	if (hypot(x-eta1o,y-eta2o) > MAXIMUM_DELTA_ETA) continue;
	if (tried.insert(*j).second) {
	  // Evaluate if have not already tried this position:
	  processSample(*j); 
	  if (j->lnProb < minProb) minProb = j->lnProb;
	  if (j->lnProb < minProb + PROB_THRESHOLD) {
	    gotOne = true;
	    current.push_back(*j);
	  }
	  //**/cerr << x << " " << y << " " << exp(j->lnProb) << " " << j->level << endl;
	}
      }
    }
    //cout << "shape2 did it while gotOne" << endl;
  } while (gotOne);

  // Clean out sub-threshold points
  for (list<Sample>::iterator i=current.begin(); i!=current.end(); ) {
    if (i->lnProb >= minProb + PROB_THRESHOLD) i=current.erase(i);
    else ++i;
  }

  double dA = sqrt(3.)/2. * pow( MINSTEP * pow(2., START_LEVEL), 2.);

  sampleDensity = MINSTEP * pow(2., START_LEVEL);

  for (int spawnLevel = START_LEVEL - 1; spawnLevel>=0; spawnLevel--) {

    /**cerr << "After spawn level " << spawnLevel+1
	     << " trials " << tried.size()
	     << " keepers " << current.size()
	     << endl;
    if (current.size()==1) {
      double x,y; current.front().getXY(x,y);
      cerr << "keeping " << x << " , " << y
	   << " prob " << current.front().lnProb
	   << endl;
    }
    */

    // Done if we already have enough points
    if (current.size() >= MINPOINTS) break;
    //**/cerr << "**At spawn level " << spawnLevel << endl;

    dA /= 4.;
    sampleDensity /= 2.;

    list<Sample> spawners = current;
    for (list<Sample>::iterator i=spawners.begin(); i!=spawners.end(); ++i) {
      // Skip this guy if he's sub-threshold
      if (i->lnProb >= minProb + PROB_THRESHOLD) continue;
      // Get its spawnees
      list<Sample> tryem = i->neighbors(spawnLevel);
      for (list<Sample>::iterator j=tryem.begin(); j!=tryem.end(); ++j) {
	double x,y;
	j->getXY(x,y);
	// Skip this one if out of bounds
	if (hypot(x-eta1o,y-eta2o) > MAXIMUM_DELTA_ETA) continue;
	// Skip this one if we already tried it
	if (!tried.insert(*j).second) continue;
	// Try this one
	processSample(*j);
	if (j->lnProb < minProb) minProb = j->lnProb;
	// Add to accepted list if above threshold
	if (j->lnProb < minProb + PROB_THRESHOLD) current.push_back(*j);
	//**/cerr << x << " " << y << " " << exp(j->lnProb) << " " << j->level << endl;
      } // done with this member's spawn
    } // Done with new spawn.

    // Clean out sub-threshold points
    for (list<Sample>::iterator i=current.begin(); i!=current.end();) {
      if (i->lnProb >= minProb + PROB_THRESHOLD) i=current.erase(i);
      else ++i;
    }
  }

  // We should now have our list of above-threshold points.  Dump 
  double sump = 0.;
  double sumx = 0.;
  double sumxx = 0.;
  double sumy = 0.;
  double sumyy = 0.;
  double sumxy = 0.;
  double sumpCentroidMismatch=0.;
  double sumpTooSmall=0.;
  double sumpTooBig=0.;
  double sumTooElliptical=0.;

#define SAMPLING

#ifdef SAMPLING
  double sumxxx = 0.;
  double sumxxxx = 0.;
  double sumyyy = 0.;
  double sumyyyy = 0.;
#endif

  for (list<Sample>::iterator i=current.begin(); i!=current.end(); ++i) {
    double prob = exp(-i->lnProb);
    double x, y;
    i->getXY(x,y);
    sump += prob;
    sumx += prob*x;
    sumy += prob*y;
    sumxx += prob*x*x;
    sumyy += prob*y*y;
    sumxy += prob*x*y;
#ifdef SAMPLING
    sumxxx += prob*x*x*x;
    sumxxxx += prob*x*x*x*x;
    sumyyy += prob*y*y*y;
    sumyyyy += prob*y*y*y*y;
#endif
    if (!centerIsFixed && (i->centroidFlag & CentroidMismatch)) sumpCentroidMismatch+=prob;
    if (!sizeIsFixed) {
      sumpTooSmall += prob * i->fractionTooSmall;
      sumpTooBig += prob * i->fractionTooBig;
    }
    //if (hypot(x,y) > WARN_ETA_THRESHOLD) sumTooElliptical += prob;
  }

  // Set centroid flag if it was set for significant fraction
  // of the likelihood
    /*/
    cerr << "fraction of likelihood which has centroid mismatch flag set: "
         << (sumpCentroidMismatch / sump) << " >? " << CENTROID_FLAG_THRESHOLD << endl;
    /*/
  if (sumpCentroidMismatch / sump > CENTROID_FLAG_THRESHOLD) {
    flags |= CentroidMismatch;
  } else {
    flags &= ~CentroidMismatch;
  }

  if (!sizeIsFixed) {
    // If too much likelihood is piled up at the bounds of
    // the sigma range, fix the galaxy sigma at the limit
    // and repeat the process
    if ( sumpTooSmall / sump > SIZE_PILEUP_THRESHOLD) {
      // set size to minimum
      flags |= SizeFixed;
      galaxyBasis.setMu(lnSigmaFloor);
      setSizing(false);
    } else if ( sumpTooBig / sump > SIZE_PILEUP_THRESHOLD) {
      // set size to maximum
      flags |= SizeFixed;
      galaxyBasis.setMu(lnSigmaCeiling);
      setSizing(false);
    }
    if (sizeIsFixed) return shape2(logLikelihood, covE);
  }

  // Set flag if not enough points over likely region
  if (current.size() < MINPOINTS){
    //cerr << "only " << current.size() << " points over likely region: undersampled!" << endl;
    flags |= UnderSampled;
  }
  else {
    //cerr << "a good " << current.size() << " points over likely region." << endl;
    flags &= ~UnderSampled;
  }

  // Set flag if too much probability is at large eta
  //if (sumTooElliptical / sump > WARN_ETA_PROBABILITY) {
  //  flags |= TooElliptical;
  //} else {
    flags &= ~TooElliptical;
  //}

  eTrialCount = tried.size();
  totalProbability = sump * dA / (2*PI);

  Shear meanS;
  double eta1 = sumx / sump;
  double eta2 = sumy / sump;
  meanS.setEta1Eta2(eta1,eta2);
  covE(0,0) = MAX(0.,sumxx/sump - eta1*eta1);
  covE(0,1) = sumxy/sump - eta1*eta2;
  covE(1,1) = MAX(0.,sumyy/sump - eta2*eta2);
  logLikelihood = minProb;

#ifdef SAMPLING
  double sx=0., sy=0., sp=0.;
  for (list<Sample>::iterator i=current.begin(); i!=current.end(); ++i) {
    double x, y;
    i->getXY(x,y);
    DVector dx(2);
    dx[0] = x-eta1; dx[1] = y-eta2;
    double prob = exp(-0.5*dx*(dx/i->covE))/sqrt(i->covE.det());
    sx += prob*dx[0];
    sy += prob*dx[1];
    sp += prob;
  }

  double skewx = (sumxxx/sump - 3*eta1*covE(0,0) - pow(eta1,3.))*pow(covE(0,0), -1.5);
  double skewy = (sumyyy/sump - 3*eta2*covE(1,1) - pow(eta2,3.))*pow(covE(1,1), -1.5);
  double kurtx = (sumxxxx/sump - 4*sumxxx*eta1/sump + 6*covE(0,0)*eta1*eta1 + 
		  3*pow(eta1, 4.)) * pow(covE(0,0),-2.) - 3.;
  double kurty = (sumyyyy/sump - 4*sumyyy*eta2/sump + 6*covE(1,1)*eta2*eta2 + 
		  3*pow(eta2, 4.)) * pow(covE(1,1),-2.) - 3.;
  /*
  cerr << "#Skew " << skewx << " " << skewy << " kurt " << kurtx << " " << kurty  
       << " debias " << sx/sp << " " << sy/sp
       << endl;

  for (int isamp=0; isamp<10; isamp++) {
    double p=u * sump;
    for (list<Sample>::iterator i=current.begin(); i!=current.end(); ++i) {
      p -= exp(-i->lnProb);
      double x, y;
      i->getXY(x,y);
      if (p<0) {
	cerr << "#Sample " << x << " " << y << endl;
	break;
      }
    }
  }
  */
#endif
    
  return meanS;
}

///////////
void
fitQuadratic(int index, DVector& a, double& upperLimit,
	     list<double>& x,
	     list<double>& y) {
  Assert(x.size()==y.size());
  if (x.size() < 2) throw FDNTError("fitQuadratic attempted with <2 points");
  if (x.size()==2) {
    // Only 2 points, do a linear fit
    if (index>0) throw FDNTError("fitQuadratic with 2 points, starting past 1st pt");
    upperLimit = x.back();
    DMatrix X(2,2, 1.);
    DVector Y(2);
    list<double>::iterator xp=x.begin();
    list<double>::iterator yp=y.begin();
    int i=0;
    while (xp!=x.end()) {
      X(i,1) = *xp;
      Y[i] = *yp;
      ++xp;
      ++yp;
    }
    Y /= X;
    a[0] = Y[0]; a[1]=Y[1]; a[2]=0.;
    return;
  }
  // Quadratic fit:
  if (index > x.size()-3) 
    FormatAndThrow<FDNTError>() << "fitQuadratic with index " << index
				<< " and size " << x.size();
  DMatrix X(3,3,1.);
  list<double>::iterator xp=x.begin();
  list<double>::iterator yp=y.begin();
  for (int i = 0; i<index; i++) {
    ++xp; ++yp;
  }
  for (int i=0; i<3; i++) {
    X(i,1) = *xp;
    X(i,2) = *xp * *xp;
    a[i] = *yp;
    ++xp;
    ++yp;
  }
  a /= X;
  if (xp==x.end()) upperLimit = x.back();
  else {
    --xp;
    upperLimit = *xp;
    --xp;
    upperLimit = 0.5*(upperLimit + *xp);
  }
  return;
}

// Probability only, includes numerical marginalization
// over ln(sigma)
template <class T>
double 
FDNT<T>::logProbability2(Shear targetS,	
			 double& fractionTooSmall,
			 double& fractionTooBig,
			 tmv::SymMatrix<double>& covE)
{
  
  // Standard deviation of ln(sigma) that is small enough to assume Gaussian 
  // integration is valid:
  const double ANALYTIC_SIZE_ERROR=0.1;

  // How far to integrate beyond above points to check
  // for substantial probability beyond bound
  const double INTEGRATION_BUFFER = 0.4;
  // Step size in ln(sigma) to calculate lnProb:
  const double LNSIG_STEP_INIT = 0.2;
  // Step size for numerical integration with interpolation
  const double LNSIG_INTEGRATION_STEP_INIT=0.01;
  // Factor to divide by for refining ln(sigma) step sizes
  const double LNSIG_DIV_FACTOR=10.;
  // Probability threshold for ln(sigma)
  const double LNSIG_PROB_THRESHOLD=-log(0.002)*2;

  const double lnSigmaFloor = log(leastPSFSigma)-TOO_SMALL;
  const double lnSigmaCeiling = nativeBasis.getMu() + TOO_BIG;

  // Continue analytic marginalization until ln(likelihood)
  // is this far away from most likely point:
  const double LNLIKELIHOOD_CUTOFF = 5.;
  // Map ln(sigma) likelihood until ends are this far above peak:
  const double LNPROB_DEPTH = 3.;

  int iE1, iE2, iSize, iX, iY, nRe, nIm, nTests;
  //cout << "logProbability2: settingIndices" << endl;
  ExposureGroup<float>::setIndices(centerIsFixed, sizeIsFixed,
				   iE1, iE2, iSize, iX, iY, nRe, nIm, nTests);

  fractionTooSmall = 0.;
  fractionTooBig = 0.;

  tmv::SymMatrix<double> covE23(nRe, 0.);
  DVector dE23(nRe, 0.);
  //cout << "logProbability2: calling marginalizedCentroid" << endl;
  marginalizedCentroid(targetS, dE23, covE23);
  //cout << "logProbability2: called marginalizedCentroid" << endl;
  
  if (iSize<0) {
    //cout << "logProbability2: // No centroid marginalization necessary:" << endl;
    covE23.saveDiv();
    double chisq = dE23 * (dE23 / covE23);
    covE = covE23;
    return 0.5 * ( chisq + log(covE23.det()) );
  }

  // If max likelihood ln(sigma) is close to test value, and
  // the Fisher matrix is narrow in ln(sigma), then use
  // analytic Gaussian marginalization

  // Calculate the mu shift needed to minimize lnProb at this shear:
  double gaussShift = -dE23[iSize] 
    + covE23.col(iSize, iE1,iE2+1) * (dE23.subVector(iE1, iE2+1) / covE23.subSymMatrix(iE1, iE2+1));
    
  //cout <<  "logProbability2: got gaussShift " << gaussShift << endl;
    
  if (false) { /**/
    //  if (abs(gaussShift) < SIZE_MISMATCH_THRESHOLD
    //      && covE23(iSize,iSize) < ANALYTIC_SIZE_ERROR*ANALYTIC_SIZE_ERROR) {
    tmv::SymMatrix<double> covE2 = covE23.subSymMatrix(iE1, iE2+1);
    DVector dE2 = dE23.subVector(iE1, iE2+1);
    covE2.saveDiv();
    double chisq = dE2 * (dE2 / covE2);
    double lnProb = 0.5 * ( chisq + log(covE2.det()) );
    // fractionTooSmall = ????
    /**/cerr << targetS << "Gaussian marginalization " << galaxyBasis.getMu() + gaussShift
	     << " covSize " << sqrt(covE23(iSize,iSize))
	     << " FishSize " << pow( covE23.inverse()(iSize,iSize), -0.5)
	     << " lnprob " << lnProb
	     << endl;
    /**/cerr << "covE23:\n" << covE23 << "dE23:\n" << dE23 << endl;
    //**??    galaxyBasis.setMu( galaxyBasis.getMu() + gaussShift );
    //**??    return lnProb;
  }

  // We will have to map out the ln(sigma) likelihood numerically
  // Make a list of calculated values
  list<double> vlnsig;
  list<double> vlnprob;

  // Put current probability in the array and set it as min
  covE23.saveDiv();
  //cout << "logProbability2: savedDiv" << endl;
  double chisq = dE23 * (dE23 / covE23);
  double lnprob = 0.5 * ( chisq + log(covE23.det()) );

  if(isnan(lnprob) || isinf(lnprob))
  {
    flags |= Singularity;
    FormatAndThrow<FDNTError>() << "Invalid probability encountered in numerical marginalization";
  }


  covE23.unsetDiv();
  //cout << "logProbability2: unsetDiv" << endl;

  /** tmv::SymMatrix<double> gcov=covE23; gcov.saveDiv();
   DVector gde = dE23;
   double glnsig = galaxyBasis.getMu() - dE23[iSize];
   **/
  double minlnprob = lnprob;
  double minlnprobsig = galaxyBasis.getMu();

  bool lnsigStepTooLarge;
  double lnsigStep = LNSIG_STEP_INIT;
  double lnsigIntegrationStep=LNSIG_INTEGRATION_STEP_INIT;

  do {
    lnsigStepTooLarge = false;
    vlnsig.clear();
    vlnprob.clear();
    // add the presumed minimum point as first element
    vlnsig.push_back(minlnprobsig);
    vlnprob.push_back(minlnprob);
  
    //cout << "logProbability: pushed back galaxy mu" << endl;
    //cerr << "lnsig step size is: " << lnsigStep << endl;

    // First move is up or down depending upon sign of DE:
    bool goUp = dE23[iSize] < 0;

    // do block for iterating around the lnprob minimum
    do {
      //cout << "logProbability2: // Calculate new point" << endl;
      double nextlnsig = goUp ? 
	vlnsig.back() + lnsigStep :
	vlnsig.front() - lnsigStep ;
      galaxyBasis.setMu(nextlnsig);
      marginalizedCentroid(targetS, dE23, covE23);
      covE23.saveDiv();
      chisq = dE23 * (dE23 / covE23);
      lnprob = 0.5 * ( chisq + log(covE23.det()) );
      covE23.unsetDiv();
      // check for lnsig step size
      if (abs(lnprob - minlnprob) > LNSIG_PROB_THRESHOLD) lnsigStepTooLarge = true; 
      // Mark new minimum
      if (lnprob < minlnprob) {
	minlnprob = lnprob;
	minlnprobsig = nextlnsig;
      }
      // And save into ordered array
      if (goUp) {
	vlnsig.push_back(nextlnsig);
	vlnprob.push_back(lnprob);
      } else {
	vlnsig.push_front(nextlnsig);
	vlnprob.push_front(lnprob);
      } 

      // Now decide which direction next trial goes,
      // or possibly finished
      
      bool downDone = vlnsig.front() < lnSigmaFloor-INTEGRATION_BUFFER 
	|| vlnprob.front() > minlnprob + LNPROB_DEPTH;
      bool upDone = vlnsig.back() >  lnSigmaCeiling + INTEGRATION_BUFFER
	|| vlnprob.back() > minlnprob + LNPROB_DEPTH;
      
      // Finished?
      if (upDone && downDone) break; 
    
      // Otherwise, go up if down is done, vice-versa
      if (upDone) {
	goUp = false;
      } else if (downDone) {
	goUp = true;
      } else {
	// Go whichever direction has most likelihood
	goUp = vlnprob.back() < vlnprob.front();
      }
    } while (true);	// ??? put some limit on iterations??

    galaxyBasis.setMu(minlnprobsig);
    marginalizedCentroid(targetS, dE23, covE23); // need dE23 if lnsigStepTooLarge

    lnsigStep /= LNSIG_DIV_FACTOR;
    lnsigIntegrationStep /= LNSIG_DIV_FACTOR;

  } while (lnsigStepTooLarge);
  

  /**   // DEBUG CODE
  static const char alphanum[] =
    "0123456789"
    "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    "abcdefghijklmnopqrstuvwxyz";
  int fnlen = 4;
  char fname[fnlen+1];
  for (int i=0; i<fnlen; ++i) {
    fname[i] = alphanum[rand() % (sizeof(alphanum)-1)];
  }
  fname[fnlen] = 0;
  ofstream ofs(fname);
  ofs << "#lnsig   lnprob" << endl;

  list<double>::iterator j=vlnprob.begin(); 
  for (list<double>::iterator i=vlnsig.begin();
       i!=vlnsig.end(); 
       ++i, ++j) {
    ofs << *i << " " << *j << endl;
  }
  ofs.close();
  **/

  // Now take the array of probabilities and integrate, fitting
  // quadratic to vicinity of each point
  double integrateMin = lnSigmaFloor - INTEGRATION_BUFFER;
  double integrateMax = lnSigmaCeiling + INTEGRATION_BUFFER;

  double sumProb = 0;
  double sumProbTooSmall = 0;
  double sumProbTooBig = 0;
  double sumProbLnsig = 0.;

  // ??? Do something if Max<=Min

  /**/
  if (false) {
    static int ctr=0;
    ctr++;
    // Find index of max-likelihood point
    int minIndex=0;
    for (list<double>::iterator yp=vlnprob.begin();
	 yp != vlnprob.end() && *yp > minlnprob+0.01;
	 ++yp, ++minIndex) ;
    if (minIndex>0) minIndex--;
    DVector a(3);	// quadratic coefficients
    double fitUpperLimit;	// Largest lnsig to use with this fit
    fitQuadratic(minIndex, a, fitUpperLimit, vlnsig, vlnprob);
    DVector alin(3);	// quadratic coefficients fit to linear sigma
    list<double> vsig;
    list<double> linprob;
    list<double>::iterator yp=vlnprob.begin();
    for (list<double>::iterator xp=vlnsig.begin();
	 xp != vlnsig.end();
	 ++xp, ++yp) {
      vsig.push_back(exp(*xp));
      linprob.push_back(*yp  + *xp);
    }
    fitQuadratic(minIndex, alin, fitUpperLimit, vsig, linprob);

    list<double>::iterator xp=vlnsig.begin();
    yp=vlnprob.begin();
    for ( ; xp!=vlnsig.end(); ++xp, ++yp)  {
      double sig=exp(*xp);
      cerr << ctr << " " << *xp << " " << *yp 
	   << " logquad " << a[0] + *xp*(a[1] + a[2]* *xp)
	   << " linquad " << (alin[0] + sig*(alin[1] + alin[2]*sig)) - *xp
	   << endl;
    }
  } /**/

  // Start by fitting quadratic to first points
  int fitIndex=0;
  DVector a(3);	// quadratic coefficients
  double fitUpperLimit;	// Largest lnsig to use with this fit
  fitQuadratic(fitIndex, a, fitUpperLimit, vlnsig, vlnprob);

  for (double lnsig = integrateMin;
       lnsig <= integrateMax;
       lnsig += lnsigIntegrationStep) {
    while (lnsig < vlnsig.back() && lnsig > fitUpperLimit) {
      fitIndex++;
      fitQuadratic(fitIndex, a, fitUpperLimit, vlnsig, vlnprob);
    }
    lnprob = lnsig*(lnsig*a[2]+a[1]) + a[0];
    // Don't go beyond cutoff for integration:
    if (lnprob > minlnprob + LNLIKELIHOOD_CUTOFF) continue;
    /** if (false)  {
      galaxyBasis.setMu(lnsig);
      marginalizedCentroid(targetS, dE23, covE23);
      covE23.saveDiv();
      chisq = dE23 * (dE23 / covE23);
      double lnprob2 = 0.5 * ( chisq + log(covE23.det()) );
      covE23.unsetDiv();
      gde[iSize] = lnsig - glnsig;
      cerr << lnsig << " " << lnprob 
	       << " gauss " << 0.5*(gde*(gde/gcov) + log(gcov.det()))
	       << " real " << lnprob2
	       << endl;
    }
    **/

    double prob = exp( minlnprob - lnprob);
    sumProb += prob;
    sumProbLnsig += prob*lnsig;
    if (lnsig < lnSigmaFloor)
      sumProbTooSmall += prob;
    if (lnsig > lnSigmaCeiling + TOO_BIG)
      sumProbTooBig += prob;
  }

  if (sumProb == 0.) {
    // only probability is outside bounds??? 
    // this should be rare/odd.
    return 999.;
  }

  fractionTooSmall = sumProbTooSmall / sumProb;
  fractionTooBig = sumProbTooBig / sumProb;
  galaxyBasis.setMu(sumProbLnsig / sumProb);
  /**cerr << targetS << " Numerical marginalization "
	   << " at mean ln(sigma) " << galaxyBasis.getMu()
	   << " lnprob " << minlnprob - log(sumProb * lnsigIntegrationStep) + 0.5*log(2*PI)
	   << endl;
  cerr << " ... minlnprob " << minlnprob
       << " integral " << sumProb * lnsigIntegrationStep
       << " a " << a[0] << " " << a[1] << " " << a[2] << endl;
	   **/
  // ???? set covE here ???
  return minlnprob - log(sumProb * lnsigIntegrationStep) + 0.5*log(2*PI);
}

template <class T>
void
FDNT<T>::wtFlux(double& f, double& varf) {
  double ft, varft;
  f = varf = 0.;
  for (int i=0; i<veg.size(); i++) {
    if (!useEG[i]) continue;
    veg[i]->wtFlux(galaxyBasis.getS(), *wg, ft, varft);
    f += ft / varft;
    varf += 1./ varft;
  }
  f /= varf;
  varf = 1/varf;
}


////////////////////////////////////////////////////
/// Iteration methods
////////////////////////////////////////////////////

#if 0
/***** ???? revive later

// Here's the routine to iterate toward zero null tests:
template <class T>
bool 
FDNT<T>::doglegNull(double etaTol) {
  DVector eta(2);
  DVector M(2);
  DMatrix dMdE(2,2);
  // Start at whatever current estimate of galaxyBasis is:
  Shear targetS = galaxyBasis.getS();

  // Initial trial - want real dMdE's here to seek null:
  targetS.getEta1Eta2(eta[0],eta[1]);
  try {
    M = tests(targetS, dMdE);
  } catch (FDNTError& e) {
    return false;
  }

  double trustRadius = INITIAL_TRUST_RADIUS;
  int istep;
  bool success = false;

#ifdef DEBUG
  cerr << "Dogleg with sigma " << exp(galaxyBasis.getMu()) 
       << " and native basis " << nativeBasis << endl;
#endif
  for (istep=0; istep<MAX_ITERATIONS; istep++) {
#ifdef DEBUG
    cerr << ".. step " << istep << " " << eta << endl;
#endif
    {
      double normeta = hypot(eta[0],eta[1]);
      if (normeta > MAXIMUM_ETA) {
	eta /= normeta;
	targetS.setEta1Eta2(eta[0],eta[1]);
	galaxyBasis.setS(targetS);
	flags |= TooElliptical;
	return false;
      }
    }
    // Quit if eta is too large:
    DVector nr(2,0.);
    DVector dE(2,0.);
    if (dMdE.det() == 0.) {
      flags |= Singularity;
      return false;
    }

    // Estimate Newton-Raphson step, use if it's short enough
    nr = - M/dMdE;
    if ( hypot(nr[0],nr[1]) < trustRadius) {
      dE = nr;
    } else {
      // Calculate steepest-descent step
      DVector grad = -dMdE.transpose()*M;
      DVector tmp = dMdE*grad;
      grad *= (grad[0]*grad[0]+grad[1]*grad[1])  /
	(tmp[0]*tmp[0]+tmp[1]*tmp[1]);
      if ( hypot(grad[0],grad[1])>trustRadius) {
	// Follow steepest descent to trust radius:
	double shrink = trustRadius / hypot(grad[0],grad[1]);
	dE = grad * shrink;
      } else {
	// Find where dogleg crosses trust radius
	DVector b=nr - grad;
	double aa = grad*grad;
	double ab = grad*b;
	double bb= b*b;
	double alpha = (-ab+sqrt(ab*ab+bb*(trustRadius*trustRadius-aa)))/bb;
	if (alpha<0. || alpha>1.) {
	  cerr << "Bad alpha; nr=" << nr
	       << " grad " << grad
	       << " trustRadius " << trustRadius
	       << " alpha " << alpha
	       << endl;
	  exit(1);
	}
	dE = grad + alpha*b;
      }
    }

    // If the recommended step is short enough, we have "succeeded"
    if (hypot(dE[0],dE[1]) < etaTol) {
	targetS.setEta1Eta2(eta[0],eta[1]);
	galaxyBasis.setS(targetS);
	return true;
    }

    DVector expecteddM = dMdE * dE;
    double expecteddMsq = -2*(M*expecteddM) - expecteddM*expecteddM;

    // Make a trial at new point
#ifdef DEBUG
    cerr << " trial dE " << dE << endl;
#endif
    targetS.setEta1Eta2(eta[0]+dE[0], eta[1]+dE[1]);
    DMatrix newdMdE(2,2);
    DVector newM = tests(targetS, newdMdE);
    double gainFactor = (M[0]*M[0]+M[1]*M[1] - newM[0]*newM[0] -newM[1]*newM[1])
      / expecteddMsq;

    // Was there improvement?
    if (gainFactor > 0.) {
      eta += dE;
      M = newM;
      dMdE = newdMdE;
    } 

    // Adjust the trustRadius
    if (gainFactor < 0.25) {
      trustRadius *= 0.5;
      if (trustRadius < MINIMUM_TRUST_RADIUS) {
	flags |= DidNotConverge;
#ifdef DEBUG
	cerr << "Under trust, gain factor " << gainFactor << endl;
#endif
	targetS.setEta1Eta2(eta[0],eta[1]);
	galaxyBasis.setS(targetS);
	return false;
      }
    } else if (gainFactor > 0.75) {
      trustRadius = MAX(trustRadius, 3.*hypot(dE[0],dE[1]));
      trustRadius = MIN(trustRadius, INITIAL_TRUST_RADIUS);
    }
  }

  // If we run to max steps, we have failed
  flags |= DidNotConverge;
  targetS.setEta1Eta2(eta[0],eta[1]);
  galaxyBasis.setS(targetS);
  return false;
}


// Line minimization class
template <class T>
class FDNTLineFunction {
private:
  FDNT<T>& fd;
  DVector E0;
  DVector dE;
  tmv::SymMatrix<double>& covE;
public:
  FDNTLineFunction(FDNT<T>& fd_, 
		   const DVector& E0_, 
		   const DVector& dE_,
		   tmv::SymMatrix<double>& covE_): fd(fd_), E0(E0_), dE(dE_),
						   covE(covE_) {
    double h=hypot(dE[0],dE[1]); 
    dE /= h;
  }
  double operator()(double dEta) {
    Shear trialS;
    trialS.setEta1Eta2(E0[0]+dE[0]*dEta, E0[1]+dE[1]*dEta);
    //** cerr << "FDNTLineFunction looking at "
    //**	     <<  E0[0]+dE[0]*dEta << " " <<  E0[1]+dE[1]*dEta << endl;
    // Use the cleanDerivatives to get better-behaved probabilities
    return fd.logProbability(trialS, covE);
  }
};

template <class T>
Shear
FDNT<T>::shape(double& logLikelihood,
	       tmv::SymMatrix<double>& covE) {
  //**cerr << "Starting with mu " << galaxyBasis.getMu() << endl;
  DVector oldE(2);
  double oldProb;
  double startEtaStep=0.01;
  double etaTolerance=0.0001;
  double probTolerance=0.001;
  // Initialize start point and direction vectors
  Shear trialS = galaxyBasis.getS();
  trialS.getEta1Eta2(oldE[0], oldE[1]);

  DVector dir1(2,0.);
  DVector dir2(2,0.);
  dir1[0] = 1.;
  dir2[1] = 1.;

  //clear some flags:
  flags &= ~(TooElliptical + DidNotConverge + CentroidMismatch + SizeMismatch);

  for (int steps=0; steps<MAX_ITERATIONS; steps++) {
    // Line minimization in successive directions
    try {
      double de1, de2;
      {
	FDNTLineFunction<T> lf(*this, oldE, dir1, covE);
	if (steps==0) oldProb = lf(0.);
#ifdef DEBUG
	cerr << "***step " << steps
	     << " start " << oldE 
	     << " prob " << oldProb
	     << endl;
#endif
	brent::Brent<FDNTLineFunction<T> > br(lf, startEtaStep);
	br.setXTolerance(etaTolerance);
	double p=br.minimum(de1);
#ifdef DEBUG
	cerr << "Move " << de1
	     << " along " << dir1[0] << " " << dir1[1]
	     << " p " << p
	     << endl;
#endif
      }
      {
	FDNTLineFunction<T> lf(*this, oldE + de1*dir1, dir2, covE);
	brent::Brent<FDNTLineFunction<T> > br(lf, startEtaStep);
	br.setXTolerance(etaTolerance);
	double p=br.minimum(de2);
#ifdef DEBUG
	cerr << "Move " << de2
	     << " along " << dir2[0] << " " << dir2[1]
	     << " p " << p
	     << endl;
#endif
      }
      // No motion: done
      if (de1==0. && de2==0.) {
	// trialS and covE were set on final evaluation.
	logLikelihood = oldProb;
	trialS.setEta1Eta2(oldE[0], oldE[1]);
	galaxyBasis.setS(trialS);
	return trialS;
      }

      // Then along the "average" direction
      DVector dirSum = de1*dir1 + de2*dir2;
      DVector newE = oldE + dirSum;
      double deTry = sqrt(dirSum*dirSum);
      dirSum /= deTry;
      double newProb;
      {
	FDNTLineFunction<T> lf(*this, newE, dirSum, covE);
	brent::Brent<FDNTLineFunction<T> > br(lf, deTry);
	br.setXTolerance(etaTolerance);
	double de;
	newProb = br.minimum(de);
	newE += de*dirSum;
#ifdef DEBUG
	cerr << "Move " << de
	     << " along " << dirSum[0] << " " << dirSum[1]
	     << " p " << newProb
	     << endl;
#endif
      }
      // Quit if we have moved outside the maximum eta:
      {
	double eta = hypot(newE[0], newE[1]);
	if (eta > MAXIMUM_ETA) {
	  trialS.setEta1Eta2(MAXIMUM_ETA*newE[0]/eta, 
			     MAXIMUM_ETA*newE[1]/eta);
	  galaxyBasis.setS(trialS);
	  flags |= TooElliptical;
	  return trialS;
	}
      }

      // Quit if change in probability is to small or motion is small
      DVector moved = newE - oldE;
      
      if (abs(oldProb - newProb) < probTolerance  ||
	  moved*moved < etaTolerance*etaTolerance) {
	logLikelihood=newProb;
	trialS.setEta1Eta2( newE[0], newE[1]);
	galaxyBasis.setS(trialS);
	return trialS;
      }

      oldE=newE;
      oldProb = newProb;

      // Use average direction instead of the one most parallel to it.
      if (abs(dirSum*dir1) > abs(dirSum*dir2)) {
	dir1 = dir2;
      }
      dir2 = dirSum;
    } catch (brent::BrentError& b) {
      flags |= DidNotConverge;
      return trialS;
    }
  }
  // Too many steps
  flags |= DidNotConverge;
  return trialS;
}

template <class T>
void
FDNT<T>::marginalizedSize(Shear targetS,
			  DVector& dE2,
			  tmv::SymMatrix<double>& covE2) {
  int iE1, iE2, iSize, iX, iY, nRe, nIm, nTests;
  ExposureGroup<float>::setIndices(centerIsFixed, sizeIsFixed,
				   iE1, iE2, iSize, iX, iY, nRe, nIm, nTests);

  const int MAXIMUM_SIZE_ITERATIONS = 4;
  // Check input vector sizes
  Assert(dE2.size()==2);
  Assert(covE2.nrows()==2);

  tmv::SymMatrix<double> covE23(nRe, 0.);
  DVector dE23(nRe, 0.);
  marginalizedCentroid(targetS, dE23, covE23);

  if (iSize<0) {
    // No centroid marginalization necessary:
    dE2 = dE23;
    covE2 = covE23;
    return;
  }

  //**vector<double> ldE;
  //**vector<double> lsig;

  for (int iterationCount=1; 
       iterationCount<MAXIMUM_SIZE_ITERATIONS;
       iterationCount++) {

    //**ldE.push_back(dE23[iSize]); lsig.push_back(galaxyBasis.getMu());

    // Is max likelihood galaxy sigma far from test value?
    if (abs(dE23[iSize]) < SIZE_MISMATCH_THRESHOLD) {
      // Happy here!
      //**cerr << "Success: "  << " shear " << targetS << endl;
      //**for (int jj=0; jj<ldE.size(); jj++)
	//** cerr << "Step " << jj << " lnSig " << lsig[jj] << " dE " << ldE[jj] << endl;
      flags &= ~SizeMismatch;
      galaxyBasis.setMu( galaxyBasis.getMu() - dE23[iSize] );
      //** cerr << "Return " << galaxyBasis.getMu() 
	//**	<< " sd(lnsig) " << sqrt(covE23(iSize,iSize)) << endl;
      covE2 = covE23.subSymMatrix(iE1, iE2+1);
      dE2 = dE23.subVector(iE1, iE2+1);
      return;
    }
    // Need to re-test with different weight size.
    double lnSizeShift = dE23[iSize];
    // Cut step down to maximum size.
    if (abs(lnSizeShift) > MAXIMUM_SIZE_SHIFT)
      lnSizeShift *= MAXIMUM_SIZE_SHIFT / abs(lnSizeShift);
    double newSize = galaxyBasis.getMu() - lnSizeShift;

    // Freeze size iterations if it has gone well below the
    // PSF size or above the native size:
    bool freezeit = false;
    if (newSize < log(MINIMUM_SIZE_FACTOR*leastPSFSigma) ) {
      // Size going too small.  Freeze at lower limit
      newSize = log(MINIMUM_SIZE_FACTOR*leastPSFSigma);
      freezeit = true;
    }
    if (newSize > log(MAXIMUM_SIZE_FACTOR) + nativeBasis.getMu()) {
      // Size going too large.  Freeze at upper limit
      newSize = log(MAXIMUM_SIZE_FACTOR) + nativeBasis.getMu();
      freezeit = true;
    }
    if (freezeit) {
      //** cerr << "Freezing size at " << newSize  << " shear " << targetS << endl;
      //** for (int jj=0; jj<ldE.size(); jj++)
	//** cerr << "Step " << jj << " lnSig " << lsig[jj] << " dE " << ldE[jj] << endl;
      /** for (double ss = log(MINIMUM_SIZE_FACTOR*leastPSFSigma);
		/**	ss < log(MAXIMUM_SIZE_FACTOR) + nativeBasis.getMu();
		/**	ss += 0.2) {
	galaxyBasis.setMu(ss);
	marginalizedCentroid(targetS, dE23, covE23);
	double chisq = dE23 * (dE23 / covE23);
	double prob = 0.5 * ( chisq + log(covE23.det()) );
	cerr << "lnsig " << ss << " dE " << dE23[iSize] << " prob " << prob 
	     << " sd(lnsig) " << sqrt(covE23(iSize,iSize)) << endl;
      }

      galaxyBasis.setMu(newSize);
      setSizing(false);
      flags |= SizeFixed;
      flags &= ~SizeMismatch;
      // Calculate at fixed size and return:
      ExposureGroup<float>::setIndices(centerIsFixed, sizeIsFixed,
				       iE1, iE2, iSize, iX, iY, nRe, nIm, nTests);
      marginalizedCentroid(targetS, dE2, covE2);
      if (true) { // ???
	double chisq = dE2 * (dE2 / covE2);
	double prob = 0.5 * ( chisq + log(covE2.det()) );
	cerr << "Fixed returns prob " << prob << endl;
      } 
      return;
    }

    // Do a new measurement with updated weight size
    galaxyBasis.setMu(newSize);
    marginalizedCentroid(targetS, dE23, covE23);
  }

  //**cerr << "Overrun iterations at " << galaxyBasis.getMu()  
	//**   << " dE " << dE23[iSize] << " shear " << targetS << endl;
  //**for (int jj=0; jj<ldE.size(); jj++)
    //** cerr << "Step " << jj << " lnSig " << lsig[jj] << " dE " << ldE[jj] << endl;

  // Get here if we've overrun our size iteration count:
  flags |= SizeMismatch;
  covE2 = covE23.subSymMatrix(iE1, iE2+1);
  dE2 = dE23.subVector(iE1, iE2+1);
  return;
}

template <class T>
double 
FDNT<T>::logProbability(Shear targetS,	
		    tmv::SymMatrix<double>& covE) {
  Assert(covE.nrows()==2);

  // Get dimensions
  int iE1, iE2, iSize, iX, iY, nRe, nIm, nTests;
  ExposureGroup<float>::setIndices(centerIsFixed, sizeIsFixed,
				   iE1, iE2, iSize, iX, iY, nRe, nIm, nTests);
  DVector dE(2);
  marginalizedSize(targetS, dE, covE);
  try {
    //  Return probability
    covE.saveDiv();
    double chisq = dE * (dE / covE);
    double prob = 0.5 * ( chisq + log(covE.det()) );
    covE.unsetDiv(); // unfreeze in case user changes covE later
    return prob;
  } catch (tmv::Singular) {
    flags |= Singularity;
    covE.unsetDiv(); // unfreeze in case user changes covE later
    throw FDNTError("Singular covE matrix in probability()");
  }
}


*****/
#endif

    

  
template class FDNT<float>;
