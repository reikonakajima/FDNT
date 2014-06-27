// $Id: ExposureGroup.cpp,v 1.23 2012/02/23 22:28:01 garyb Exp $

// Code for the ExposureGroup classes of FDNT

#include "FDNT.h"
#include "fft.h"
#include "FitFlags.h"

using namespace laguerre;
using namespace fft;

// For taking all derivatives:
const double DETA=0.01;
const double DLNSHRINK=0.01;

// For storing 1/MTF and building T0:
const double MTF_FLOOR=0.001;
const double T0_CUTOFF=0.005;  // Kill filter when PSF correction is large
const double ILL_CONDITION=1e6;  // Condition number for b covariance that sets flag

// This defines the grid used for GL k-space integrations:
static const double dkGL=0.1;
static const double kmaxGL=5.;

// Size of box used for FFT:
const double FFT_BOX_SIGMA = 4.; // Default # object sigma native box extends?
const double MINIMUM_FFT_SIGMA = 3.0; // Min # object sigma native box includes?
const int MINIMUM_FFT_SIZE=16;

CircularizedSB::CircularizedSB(const sbp::SBProfile& sb): tab(GTable<>::spline),
							  maxksq(0.) {
  // tabulate the inverse-square mean of the SBProfile:
  const int NPTS=16;	// How many points for azimuthal averaging
  vector<Position<double> > cs(NPTS);
  for (int i=0; i<NPTS; i++)
    cs[i] = Position<double>(cos(i*(PI/NPTS)), sin(i*(PI/NPTS)));
  for (double k=0; k<2.*sb.maxK(); k += sb.stepK()/2.) {
    double sum=0.;
    for (int i=0; i<NPTS; i++) {
      double nn = norm(sb.kValue(cs[i] * k));
      sum += 1./MAX(nn, MTF_FLOOR*MTF_FLOOR);
      // Should really cut off when any part of annulus gets below floor ????
    }
    double t0 = 1./(sum/NPTS);
    if (t0 < T0_CUTOFF*T0_CUTOFF) {
      maxksq = k*k;
      tab.addEntry(k*k,0.);
      break;
    } else {
      tab.addEntry(k*k, t0);
    }
  }
  if (maxksq==0.) maxksq = tab.argMax();
};

template <class T>
ExposureGroup<T>::ExposureGroup(PSFInformation& psf_): 
  arrayPtr(0),
  psf(psf_),
  ktruex(0), ktruey(0),
  deconvRe(0), deconvIm(0),
  dRedxy(0), dImdxy(0),
  flags(0)
{}

template <class T>
ExposureGroup<T>::~ExposureGroup() {
  if (arrayPtr) delete arrayPtr;
  if (ktruex) delete ktruex;
  if (ktruey) delete ktruey;
  if (deconvRe) delete deconvRe;
  if (deconvIm) delete deconvIm;
  if (dRedxy) delete dRedxy;
  if (dImdxy) delete dImdxy;
}


struct laguerre::TestArrays {
public:
  const int nk;
  const int ntestsRe;
  const int ntestsIm;
  DMatrix tRe;
  DMatrix tIm;
  DMatrix tRedE;
  DMatrix tImdE;
  tmv::DiagMatrix<double> T0sq;
  tmv::DiagMatrix<double> wg;
  tmv::DiagMatrix<double> T0sqDeriv;
  tmv::DiagMatrix<double> wgDeriv;
  DMatrix kxxyy;
  DMatrix kxy;
  DVector ksq;
  bool valid(int nk_, int ntestsRe_, int ntestsIm_) {
    return nk==nk_ && ntestsRe==ntestsRe_ && ntestsIm==ntestsIm_;
  }
  TestArrays(int nk_, int ntestsRe_, int ntestsIm_):
    nk(nk_), ntestsRe(ntestsRe_), ntestsIm(ntestsIm_),
    tRe(ntestsRe, nk, 0.),
    tIm(ntestsIm, nk, 0.),
    tRedE(ntestsRe, nk, 0.),
    tImdE(ntestsIm, nk, 0.),
    T0sq(nk, 0.),
    wg(nk, 0.),
    T0sqDeriv(nk, 0.),
    wgDeriv(nk, 0.),
    kxxyy(2, nk, 0.),
    kxy(2, nk, 0.),
    ksq(nk, 0.)
  {
    //cout << "generating TestArrays with nk=" << nk_ << ", ntestsRe=" << ntestsRe_ << ", ntestsIm=" << ntestsIm_ << endl;
  }
};


template <class T>
void
ExposureGroup<T>::setIndices(bool centerIsFixed,
			     bool sizeIsFixed,
			     int& iE1, int& iE2,
			     int& iSize,
			     int& iX, int& iY,
			     int& nRe, int& nIm,
			     int& nTests) {
  // Set up indices for tests
  nTests=2;  // Always test ellipticity
  iE1 = 0;
  iE2 = 1;

  // Size test:
  iSize = sizeIsFixed ? -1 : nTests;
  if (!sizeIsFixed) nTests++;
  
  nRe = nTests;

  // Centroid tests are imaginary:
  nIm = 0;
  iX = centerIsFixed ? -1 : nTests;
  iY = centerIsFixed ? -1 : nTests+1;
  if (!centerIsFixed) nTests+=2;
 
  nIm = nTests - nRe;
}

template <class T>
DVector
ExposureGroup<T>::tests(Shear targetS,
			const CircleK& wg,
			DMatrix* dMdE,
			DVector* dMdlns,
			tmv::SymMatrix<double>* covM,
			bool centerIsFixed,
			bool sizeIsFixed,
			bool exactdMdE) {

  // *** Note that shrinkFactor and sigma_g are assumed
  // *** already set inside T0sq and wg classes.
  const int nk = ktruex->size();

  // Get indices
  int iE1,iE2,iSize,iX,iY,nRe,nIm,nTests;
  setIndices(centerIsFixed, sizeIsFixed,
	     iE1,iE2,iSize,iX,iY,nRe,nIm,nTests);
  // These are indices of x/y tests in the imaginary-valued test array
  const int iXIm = iX - nRe;
  const int iYIm = iY - nRe;
  DVector outM(nTests, 0.);

  // Get a bunch of workspace arrays if not already given them:
  if (!arrayPtr) arrayPtr = new TestArrays(nk, nRe, nIm);
  if (!arrayPtr->valid(nk, nRe, nIm)) {
    delete arrayPtr;
    arrayPtr = new TestArrays(nk, nRe, nIm);
  }

  // Some references for notational convenience:
  DMatrix& tRe = arrayPtr->tRe;
  DMatrix& tIm = arrayPtr->tIm;
  DMatrix& tRedE = arrayPtr->tRedE;
  DMatrix& tImdE = arrayPtr->tImdE;
  tmv::DiagMatrix<double>& T0sq= arrayPtr->T0sq;
  tmv::DiagMatrix<double>& vwg = arrayPtr->wg;
  tmv::DiagMatrix<double>& T0sqDeriv = arrayPtr->T0sqDeriv;
  tmv::DiagMatrix<double> wgDeriv = arrayPtr->wgDeriv;
  DVector& ksq = arrayPtr->ksq;
  DMatrix& kxxyy = arrayPtr->kxxyy;
  DMatrix& kxy = arrayPtr->kxy;

  // If we're taking derivatives, check:
  if (dMdE) {
    Assert(dMdE->ncols()==nTests && dMdE->nrows()==nTests);
    dMdE->setZero();
  }
  // If covariance matrix info requested, test size:
  if (covM)
    Assert(covM->ncols()==nTests);

  // If shrink factor deriv requested, test size:
  if (dMdlns) 
    Assert(dMdlns->size()==nTests);
  
  //cout << "tests is still alive" << endl;

  // And weight-function integrals to accumulate
  double sum_wT0=0.;
  double sum_wT0_ksq=0.;
  double sum_wDerivT0=0.;
  double sum_wDerivT0_ksq=0.;
  double sum_wT0Deriv=0.;
  double sum_wT0Deriv_ksq=0.;

  for (int i=0; i<nk; i++) {
    Position<double> ktrue((*ktruex)[i], (*ktruey)[i]);
    Position<double> ktarget = targetS.fwd(ktrue);

    double kcircsq = ktarget.x*ktarget.x+ktarget.y*ktarget.y;
    ksq[i] = kcircsq;
    kxxyy(0, i) = (ktarget.x*ktarget.x-ktarget.y*ktarget.y);
    kxxyy(1, i) = 2.*ktarget.x*ktarget.y;
    kxy(0, i) = ktarget.x;
    kxy(1, i) = ktarget.y;

    double wg_ = wg(kcircsq);
    double T0sq_ = (*(psf.T0sq))(kcircsq);
    vwg(i) = wg_;
    T0sq(i) = T0sq_;

    double T0d;
    if (dMdlns) {
      T0d = psf.T0sq->dFdlnScale(kcircsq);
      T0sqDeriv(i) = T0d;
    }
    if (iSize>=0) {
      sum_wT0 += wg_ * T0sq_ * (i==dcIndex ? 0.5 : 1.);
      sum_wT0_ksq += wg_ * T0sq_ * kcircsq;
      if (dMdE) {
	double wd = wg.dFdlnScale(kcircsq);
	wgDeriv(i) = wd;
	sum_wDerivT0 += wd * T0sq_  * (i==dcIndex ? 0.5 : 1.) ;
	sum_wDerivT0_ksq += wd * T0sq_ * kcircsq;
      }
      if (dMdlns) {
	sum_wT0Deriv += wg_ * T0d  * (i==dcIndex ? 0.5 : 1.);
	sum_wT0Deriv_ksq += wg_ * T0d * kcircsq;
      }
    }
  }



  //cout << "// Build the test matrix" << endl;
  // ??? could save a little by caching vwg * T0sq
  tRe.subMatrix(iE1, iE2+1, 0,tRe.ncols()) = kxxyy * T0sq * vwg;
  double sum_IT0;
  double sum_IT0_ksq;
  if (iSize>=0) {
    tRe.row(iSize) = T0sq * ksq * sum_wT0 - T0sq.diag() * sum_wT0_ksq;
    // This is duplicate calculation, but we'll need these for derivatives:
    if (dMdE || dMdlns) {
      sum_IT0 = T0sq.diag() * *deconvRe;
      sum_IT0_ksq = ksq * T0sq * *deconvRe;
    }
  }
  if (iX>=0) {
    tIm.subMatrix(iXIm, iYIm+1, 0,tIm.ncols()) = kxy * T0sq * vwg;
  }

  //cout << "// Calculate the values of tests" << endl;
  if (nRe>0) outM.subVector(0, nRe) = tRe * (*deconvRe);
  if (nIm>0) outM.subVector(nRe, nTests) = - tIm * (*deconvIm);
  // And covariance matrix if desired
  if (covM) *covM = covarianceOf(tRe, tIm);
  // Calculate derivatives w.r.t. centroid if needed
  if (iX>=0 && dMdE) {
    // ??? need to swap the real and imaginary parts here ???
    if (nRe>0) dMdE->subMatrix(0, nRe, iX, iY+1) = tRe * dRedxy->transpose();
    if (nIm>0) dMdE->subMatrix(nRe, nTests, iX, iY+1) = - tIm * dImdxy->transpose();
  }

  //cout << "// Calculate derivatives w.r.t. weight sigma if needed" << endl;
  if (iSize>=0 && dMdE) {
    DVector tmp = kxxyy * wgDeriv * T0sq * (*deconvRe);
    (*dMdE)(iE1, iSize) = tmp[0];
    (*dMdE)(iE2, iSize) = tmp[1];
    if (iX>=0) {
      tmp = - kxy * wgDeriv * T0sq * (*deconvIm);
      (*dMdE)(iX, iSize) = tmp[0];
      (*dMdE)(iY, iSize) = tmp[1];
    }
    (*dMdE)(iSize, iSize) = sum_IT0_ksq * sum_wDerivT0 
      - sum_IT0 * sum_wDerivT0_ksq;
  }

  // Calculate derivatives w.r.t. shrink factor if needed
  if (dMdlns) {
    DVector tmp = kxxyy * vwg * T0sqDeriv * (*deconvRe);
    (*dMdlns)[iE1] = tmp[0];
    (*dMdlns)[iE2] = tmp[1];
    if (iX>=0) {
      tmp = - kxy * vwg * T0sqDeriv * (*deconvIm);
      (*dMdlns)[iX] = tmp[0];
      (*dMdlns)[iY] = tmp[1];
    }
    if (iSize>=0) {
      DVector tmp2 = T0sqDeriv * (*deconvRe);
      (*dMdlns)[iSize] = sum_IT0_ksq * sum_wT0Deriv
	- sum_IT0 * sum_wT0Deriv_ksq
	+ (tmp2 * ksq) * sum_wT0
	- tmp2.sumElements() * sum_wT0_ksq;
    }
  }

  //cout << "// Derivs wrt ellipticity if needed" << endl;
  if (dMdE) {
    if (!exactdMdE) throw FDNTError("only implemented exactdMdE right now");

    // Make differential shears
    double eta1,eta2;
    targetS.getEta1Eta2(eta1,eta2);
    for (int ie=1; ie<=2; ie++) {
      int iEthis;
      // loop over two directions of shear
      Shear Sde;
      if (ie==1) {
	Sde.setEta1Eta2(eta1+DETA, eta2);
	iEthis = iE1;
      } else {
	Sde.setEta1Eta2(eta1, eta2+DETA);
	iEthis = iE2;
      }
      for (int i=0; i<nk; i++) {
	Position<double> ktrue((*ktruex)[i], (*ktruey)[i]);
	Position<double> ktarget = Sde.fwd(ktrue);
	
	double kcircsq = ktarget.x*ktarget.x+ktarget.y*ktarget.y;
	double wg_ = wg(kcircsq);
	double T0sq_ = (*(psf.T0sq))(kcircsq);
	
	tRedE(iE1,i) = wg_ * T0sq_ * (ktarget.x*ktarget.x-ktarget.y*ktarget.y);
	tRedE(iE2,i) = wg_ * T0sq_ * 2.*ktarget.x*ktarget.y;
	if (iSize>=0) {
	  tRedE(iSize, i) = T0sq_ * (kcircsq * sum_wT0 - sum_wT0_ksq);
	}
	if (iX>=0) {
	  tImdE(iXIm,i) = wg_ * T0sq_ * ktarget.x;
	  tImdE(iYIm,i) = wg_ * T0sq_ * ktarget.y;
	}
      }

      // Get differentials.  Note that subtracting after integration
      // would be faster, but this might have less rounding error:
      if (nRe>0) dMdE->col(iEthis,0, nRe) = ( (tRedE-tRe) * (*deconvRe) )/DETA;
      if (nIm>0) dMdE->col(iEthis,nRe, nTests) = ( (tImdE-tIm) * (*deconvIm) )/(-DETA);
    } // end E1/E2 loop
  } // end if(dMdE) for eta derivatives

  //**/cerr << "In tests: outM\n" << outM << endl;
  //**/if (covM) cerr << " covM\n" << *covM << endl;
  //**/if (dMdE) cerr << " dMdE\n" << *dMdE << endl;
  //cout << "tests are done" << endl;
  return outM;
}

template <class T>
void
ExposureGroup<T>::createdXdY() {
  Assert(ktruex);
  const int nk = ktruex->size();
  if (dRedxy) delete dRedxy;
  if (dImdxy) delete dImdxy;
  dRedxy = new DMatrix(2,nk);
  dImdxy = new DMatrix(2,nk);

  for (int i=0; i<nk; i++) {
    DComplex z( (*deconvRe)[i], (*deconvIm)[i] );
    DComplex dzdx = z * DComplex(0, (*ktruex)[i]);
    DComplex dzdy = z * DComplex(0, (*ktruey)[i]);
    (*dRedxy)(0,i) = dzdx.real();
    (*dRedxy)(1,i) = dzdy.real();
    (*dImdxy)(0,i) = dzdx.imag();
    (*dImdxy)(1,i) = dzdy.imag();
  }
}


template <class T>
void
ExposureGroup<T>::wtFlux(Shear targetS,
			 const CircleK& wg,
			 double& f, double& varf) {
  // *** Note that shrinkFactor and sigma_g are assumed
  // *** already set inside T0sq and wg classes.
  const int nk = ktruex->size();

  // Make two matrices to hold data
  DMatrix testRe(1,nk,0.);
  DMatrix testIm(0,nk,0.);	// TMV allows 0-sized matrices

  for (int i=0; i<nk; i++) {
    Position<double> ktrue((*ktruex)[i], (*ktruey)[i]);
    Position<double> ktarget = targetS.fwd(ktrue);

    double kcircsq = ktarget.x*ktarget.x+ktarget.y*ktarget.y;
    testRe(0,i) = wg(kcircsq) * (*(psf.T0sq))(kcircsq);
  }
  f = (testRe * (*deconvRe))[0];
  tmv::SymMatrix<double> mvar = covarianceOf(testRe, testIm);
  varf = mvar(0,0);
}

// Draw the real-space version of filter
template <class T>
Image<>
ExposureGroup<T>::drawFilter( Shear targetS,
			      const CircleK& wg,
			      Ellipse worldBasis,
			      bool ise1) {
  FitExposure<T> fe=felist.front();
  double xw0 = worldBasis.getX0().x;
  double yw0 = worldBasis.getX0().y;
  CrudeMap map(fe.xworld, fe.yworld, xw0, yw0);
  DMatrix dPdW = map.getdPdW();
  double xp0, yp0;
  map.worldToPix(xw0, yw0, xp0, yp0);
  int ipx0 = static_cast<int> (floor(xp0+0.5));
  int ipy0 = static_cast<int> (floor(yp0+0.5));
  double dxpix = xp0 - ipx0;
  double dypix = yp0 - ipy0;

  // Get the size of the box that will hold the object:
  Ellipse pixBasis = map.worldToPix(worldBasis);
  Bounds<double> db = pixBasis.range();
  double boxSize = MAX(db.getXMax()-db.getXMin(),
		       db.getYMax()-db.getYMin())
    * 2 * FFT_BOX_SIGMA;
  // Make FFT either 2^n or 3x2^n
  int N;
  {
    double log2n = log(2.)*ceil(log(boxSize)/log(2.));
    double log2n3 = log(3.) 
      + log(2.)*ceil((log(boxSize)-log(3.))/log(2.));
    log2n3 = MAX(log2n3, log(6.));	// must be even number
    N = static_cast<int> (ceil(exp(MIN(log2n, log2n3))-0.01));
    N = MAX(N, MINIMUM_FFT_SIZE);
  }
  
  // The desired range of image data:
  Bounds<int> bfft(-N/2, N/2-1,
		   -N/2, N/2-1);

  // Draw k space filter
  double dk = 2*PI / N;
  KTable kt(N,dk);

  double sf = FDNT<T>::shrinkFactor(targetS);
  psf.T0sq->setScaleFactor(sf);

  for (int iy=0; iy<=N/2; iy++) {
    double kpixy = iy * dk;
    for (int ix=0; ix<=N/2; ix++) {
      double kpixx = ix*dk;
      double kwx = dPdW(0,0)*kpixx + dPdW(1,0)*kpixy;
      double kwy = dPdW(0,1)*kpixx + dPdW(1,1)*kpixy;
      Position<double> ktrue(kwx,kwy);
      Position<double> ktarget = targetS.fwd(ktrue);
      double kcircsq = ktarget.x*ktarget.x+ktarget.y*ktarget.y;
      double f1 = wg(kcircsq) * (*(psf.T0sq))(kcircsq);
      f1 *= ise1 ? ktarget.x*ktarget.x-ktarget.y*ktarget.y :
	2*ktarget.x*ktarget.y;
      DComplex mtf = ExposureGroup<T>::psf.T.kValue(ktrue);
      DComplex invmtf = (norm(mtf) < MTF_FLOOR*MTF_FLOOR) ? 0. : 1./mtf;
      kt.kSet(ix,iy,f1*invmtf);
    }
  }

  // Now the negative ky's - do not use ix=0 or N/2
  for (int iy=-N/2+1; iy<0; iy++) {
    double kpixy = iy * dk;
    for (int ix=1; ix<N/2; ix++) {
      double kpixx = ix*dk;
      double kwx = dPdW(0,0)*kpixx + dPdW(1,0)*kpixy;
      double kwy = dPdW(0,1)*kpixx + dPdW(1,1)*kpixy;
      Position<double> ktrue(kwx,kwy);
      Position<double> ktarget = targetS.fwd(ktrue);
      double kcircsq = ktarget.x*ktarget.x+ktarget.y*ktarget.y;
      double f1 = wg(kcircsq) * (*(psf.T0sq))(kcircsq);
      f1 *= ise1 ? ktarget.x*ktarget.x-ktarget.y*ktarget.y :
	2*ktarget.x*ktarget.y;
      DComplex mtf = ExposureGroup<T>::psf.T.kValue(ktrue);
      DComplex invmtf = (norm(mtf) < MTF_FLOOR*MTF_FLOOR) ? 0. : 1./mtf;
      kt.kSet(ix,iy,f1*invmtf);
    }
  }
  XTable* kx = kt.transform();

  Image<> im(bfft);
  for (int ix=bfft.getXMin(); ix<=bfft.getXMax(); ix++)
    for (int iy=bfft.getYMin(); iy<=bfft.getYMax(); iy++)
      im(ix,iy) = kx->xval(ix,iy);

  cerr << "write image with bounds " << bfft << " dx " << kx->getDx() << endl;
  delete kx;
  return im;
}

///////////////////////////////////////////////////////////////
// Derived class that does FT with FFT from the observed data
///////////////////////////////////////////////////////////////

template <class T>
ExposureGroupFT<T>::ExposureGroupFT(PSFInformation& psf_):
  ExposureGroup<T>(psf_),
  varRe(0)
{}

template <class T>
ExposureGroupFT<T>::~ExposureGroupFT() {
  if (varRe) delete varRe;
}

template <class T>
double
ExposureGroupFT<T>::prepare(Ellipse startBasis) {
  if (ExposureGroup<T>::felist.size()!=1)
    throw FDNTError("Can only do FT method on single-image ExposureGroups");
  FitExposure<T> fe=ExposureGroup<T>::felist.front();
  // Start by measuring the size of this object on this image
  // Measure GLexpansion of this galaxy

  const int initialOrder=4;	// Order of GL fit used to set basis
  GLSimple<T> gl(ExposureGroup<T>::felist, startBasis, initialOrder);
  gl.setCentering(false);	// Hold to common center of all images
  ExposureGroup<T>::phaseCenter = startBasis.getX0();  // Save the origin used for FT.

  if (!gl.solve()) {
    // GL failure invalidates this image
    ExposureGroup<T>::flags = gl.getFlags();
    ExposureGroup<T>::flags |= GLFailure;
    return 0.;
  }

  Ellipse nativeBasis = gl.getBasis();

  // Return the significance of GL flux
  double f, varf;
  gl.b00(f, varf);

  double xw0 = nativeBasis.getX0().x;
  double yw0 = nativeBasis.getX0().y;
  CrudeMap map(fe.xworld, fe.yworld, xw0, yw0);
  DMatrix dPdW = map.getdPdW();
  double xp0, yp0;
  map.worldToPix(xw0, yw0, xp0, yp0);
  int ipx0 = static_cast<int> (floor(xp0+0.5));
  int ipy0 = static_cast<int> (floor(yp0+0.5));
  double dxpix = xp0 - ipx0;
  double dypix = yp0 - ipy0;

  // Get the size of the box that will hold the object:
  Ellipse pixBasis = map.worldToPix(nativeBasis);
  Bounds<double> db = pixBasis.range();
  double boxSize = MAX(db.getXMax()-db.getXMin(),
		       db.getYMax()-db.getYMin())
    * 2 * FFT_BOX_SIGMA;
  //cout << "boxSize=" << boxSize << endl;
  // Make FFT either 2^n or 3x2^n
  int N;
  {
    double log2n = log(2.)*ceil(log(boxSize)/log(2.));
    double log2n3 = log(3.) 
      + log(2.)*ceil((log(boxSize)-log(3.))/log(2.));
    log2n3 = MAX(log2n3, log(6.));	// must be even number
    N = static_cast<int> (ceil(exp(MIN(log2n, log2n3))-0.01));
    N = MAX(N, MINIMUM_FFT_SIZE);
  }
  
  //cout << "N=" << N << endl;
  
  // The desired range of image data:
  Bounds<int> bfft(-N/2+ipx0, N/2-1+ipx0,
		   -N/2+ipy0, N/2-1+ipy0);
  // Size of box containing data:
  int Nfits=N;
  //**/cerr << "N " << N << " pixBasis " << pixBasis << " want " << bfft << endl;
  // If we don't have all the desired data, see if biggest imcluded
  // square around the image is good enough
  if (!fe.sci.getBounds().includes(bfft)) {
    ExposureGroup<T>::flags |= Edge;
    // Check that we have the data
      
    Nfits = 2*(fe.sci.getBounds().getXMax()-ipx0+1);
    Nfits = MIN(Nfits, 2*(ipx0-fe.sci.getBounds().getXMin()));
    Nfits = MIN(Nfits, 2*(fe.sci.getBounds().getYMax()-ipy0+1));
    Nfits = MIN(Nfits, 2*(ipy0-fe.sci.getBounds().getYMin()));

    if (Nfits < MINIMUM_FFT_SIZE) {
      //cout << "// Not enough data included" << endl;
      ExposureGroup<T>::flags |= OutOfBounds;
      return 0.;
    }

    Bounds<int> btry(-Nfits/2+ipx0, Nfits/2-1+ipx0,
		     -Nfits/2+ipy0, Nfits/2-1+ipy0);
    if (!fe.sci.getBounds().includes(btry)) {
      //cout << "// If still a problem, we have failed (should not happen)" << endl;
      ExposureGroup<T>::flags |= OutOfBounds;
      
      return 0.;
    }
    // Border of new box must include at least some minimum size
    bfft = btry;
    for (int ix=bfft.getXMin(); ix<=bfft.getXMax(); ix++) {
      Position<double> xy1 = pixBasis.inv(Position<double>(ix, bfft.getYMin()));
      Position<double> xy2 = pixBasis.inv(Position<double>(ix, bfft.getYMax()));
      if (hypot( xy1.x, xy1.y) < MINIMUM_FFT_SIGMA
	  || hypot( xy2.x, xy2.y) < MINIMUM_FFT_SIGMA) {
	//cout << "// hypot < MINIMUM_FFT_SIGMA in x" << endl;
	ExposureGroup<T>::flags |= OutOfBounds;
	return 0.;
      }
    }
    for (int iy=bfft.getYMin(); iy<=bfft.getYMax(); iy++) {
      Position<double> xy1 = pixBasis.inv(Position<double>(bfft.getXMin(),iy));
      Position<double> xy2 = pixBasis.inv(Position<double>(bfft.getXMax(),iy));
      if (hypot( xy1.x, xy1.y) < MINIMUM_FFT_SIGMA
	  || hypot( xy2.x, xy2.y) < MINIMUM_FFT_SIGMA) {
	//cout << "// hypot < MINIMUM_FFT_SIGMA in y" << endl;
	ExposureGroup<T>::flags |= OutOfBounds;
	return 0.;
      }
    }
  }


  // Collect data in pixel box, center is close to 0,0 :
  XTable xt(N, 1.,0.);
  double sumvar = 0.;
  for (int iy=-Nfits/2; iy<Nfits/2; iy++) {
    for (int ix=-Nfits/2; ix<Nfits/2; ix++) {
      xt.xSet(ix, iy, 
	      (fe.sci( ix+ipx0, iy+ipy0 ) - fe.sky) * fe.sbScaleFactor);
      double wt = fe.wt( ix+ipx0, iy+ipy0 );
      if (wt<=0) {
	//cout << "// Can't do FT with invalid data in the field" << endl;
	ExposureGroup<T>::flags |= OutOfBounds;
	return 0.;
      }
      sumvar += fe.sbScaleFactor * fe.sbScaleFactor / wt;
    }
  }

  // Do FT: 
  KTable* kt = xt.transform();
  const double dk = kt->getDk();
  // Shift FT to place origin at nativeBasis:
  kt->translate(-dxpix, -dypix); 
  double det = dPdW.det();
  ftnorm = 2. * dk*dk / det;

  // Noise associated with each real or imag DOF
  double componentVariance = 2. * sumvar * pow(dk,4.) * pow(det,-2.);

  // Allocate and fill all arrays:
  const int nk=N*(N/2+1);
  if (ktruex) delete ktruex;
  if (ktruey) delete ktruey;
  if (deconvRe) delete deconvRe;
  if (deconvIm) delete deconvIm;
  if (varRe) delete varRe;

  // Make temporary vectors to hold useful points
  vector<double> tkx; tkx.reserve(nk);
  vector<double> tky; tky.reserve(nk);
  vector<double> ftre; ftre.reserve(nk);
  vector<double> ftim; ftim.reserve(nk);
  vector<double> var; var.reserve(nk);

  // Fill ky>=0:
  //for (int iy=0; iy<=N/2; iy++) {
  for (int iy=0; iy<N/2; iy++) {
    double kpixy = iy * dk;
    //for (int ix=0; ix<=N/2; ix++) {
    for (int ix=0; ix<N/2; ix++) {
      double kpixx = ix*dk;
      double kwx = dPdW(0,0)*kpixx + dPdW(1,0)*kpixy;
      double kwy = dPdW(0,1)*kpixx + dPdW(1,1)*kpixy;

      DComplex mtf = ExposureGroup<T>::psf.T.kValue(Position<double> (kwx,kwy));
      // Place a floor on mtf to avoid division by zero.  Throw away
      // points outside this region and assume that T0sq will always
      // kill their influence.
      if (norm(mtf) < MTF_FLOOR*MTF_FLOOR) continue;

      tkx.push_back(kwx);
      tky.push_back(kwy);
      DComplex z = kt->kval(ix, iy) * ftnorm / mtf;
      // Mark origin as self-conjugate and reduce its
      // amplitude and variance to yield proper integrations over k
      if (iy==0 && ix==0) {
	dcIndex = tkx.size()-1;
	z *= 0.5;
	var.push_back(componentVariance * 0.5 / norm(mtf));
      }	else {
	var.push_back(componentVariance / norm(mtf));
      }
      ftre.push_back(z.real());
      ftim.push_back(z.imag());
    }
  }
  
  // Now the negative ky's - do not use ix=0 or N/2
  for (int iy=-N/2+1; iy<0; iy++) {
    double kpixy = iy * dk;
    for (int ix=1; ix<N/2; ix++) {
      double kpixx = ix*dk;
      double kwx = dPdW(0,0)*kpixx + dPdW(1,0)*kpixy;
      double kwy = dPdW(0,1)*kpixx + dPdW(1,1)*kpixy;
      DComplex mtf = ExposureGroup<T>::psf.T.kValue(Position<double> (kwx,kwy));
      // Place a floor on mtf to avoid division by zero.  Throw away
      // points outside this region and assume that T0sq will always
      // kill their influence.
      if (norm(mtf) < MTF_FLOOR*MTF_FLOOR) 
      {
	//cout << "MTF too small at (ix,iy)=(" << ix << "," << iy << "); discarding" << endl;
	continue;
      }

      tkx.push_back(kwx);
      tky.push_back(kwy);
      DComplex z = kt->kval(ix, iy) * ftnorm / mtf;
      var.push_back(componentVariance / norm(mtf));
      ftre.push_back(z.real());
      ftim.push_back(z.imag());
    }
  }

  delete kt; // Done with this.

  ktruex = new DVector(tkx.size());
  std::copy(tkx.begin(), tkx.end(), ktruex->begin());
  ktruey = new DVector(tky.size());
  std::copy(tky.begin(), tky.end(), ktruey->begin());
  deconvRe = new DVector(ftre.size());
  std::copy(ftre.begin(), ftre.end(), deconvRe->begin());
  deconvIm = new DVector(ftim.size());
  std::copy(ftim.begin(), ftim.end(), deconvIm->begin());
  varRe = new tmv::DiagMatrix<double>(var.size());
  std::copy(var.begin(), var.end(), varRe->begin());

  ExposureGroup<T>::createdXdY();

  /**for (int i=0; i<ktruex->size(); i++) 
    cerr << i
	 << " " << (*ktruex)[i]
	 << " " << (*ktruey)[i]
	 << " " << (*deconvRe)[i]
	 << " " << (*deconvIm)[i]
	 << " " << (*varRe)(i)
	 << endl;
  **/
  
  return f/sqrt(varf);
}

template <class T>
tmv::SymMatrix<double> 
ExposureGroupFT<T>::covarianceOf(const DMatrix& mRe, 
				 const DMatrix& mIm) {
  const int nRe = mRe.nrows();
  const int nIm = mIm.nrows();
  const int nTests = nRe + nIm;
  tmv::SymMatrix<double> cov(nTests, 0.);
  if (nRe>0) {
    //cout << "covarianceOf: nRe>0" << endl;
    //cout << mRe << endl;
    //cout << mIm << endl;
    // Use special TMV routine to build matrix product we know is symmetric
    DMatrix tmp = mRe * (*varRe);
    //cout << "tmp: " << tmp << endl;
    tmv::SymMultMM<false>(1.,
			  tmp.view(),
			  mRe.transpose(),
			  cov.subSymMatrix(0,nRe));
  }
  //cout << "covarianceOf: this did work" << endl;
  if (nIm>0) {
    //cout << "covarianceOf: nIm>0" << endl;
    // Temporary vector with zero variance in the DC part
    tmv::DiagMatrix<double> vtmp(*varRe);
    vtmp(dcIndex) = 0.;
    DMatrix tmp = mIm * vtmp;
    tmv::SymMultMM<false>(1.,
			  tmp.view(),
			  mIm.transpose(),
			  cov.subSymMatrix(nRe,nTests));
  }
  //cout << "covarianceOf done " << nRe << " " << nIm << endl;
  return cov;
}

template <class T>
void
ExposureGroupFT<T>::setPhaseCenter(double x0, double y0) {
  // skip if we're already here:
  if (x0==ExposureGroup<T>::phaseCenter.x
      && y0==ExposureGroup<T>::phaseCenter.y) return;

  double dx = x0 - ExposureGroup<T>::phaseCenter.x;
  double dy = y0 - ExposureGroup<T>::phaseCenter.y;
  ExposureGroup<T>::phaseCenter.x = x0;
  ExposureGroup<T>::phaseCenter.y = y0;
  if (! (deconvRe || deconvIm || ExposureGroup<T>::dRedxy || ExposureGroup<T>::dImdxy ) )
    return;	// Nothing to change!
  // Build cos & sin vectors as diag matrices
  Assert(ktruex);
  const int nk=ktruex->size();
  tmv::DiagMatrix<double> C(nk, 0.);
  tmv::DiagMatrix<double> S(nk, 0.);
  DVector phase = dx * (*ktruex) + dy * (*ktruey);
  for (int i=0; i<nk; i++) {
    C(i) = cos(phase[i]);
    S(i) = sin(phase[i]);
  }
  if (deconvRe) {
    // Phase shift the Fourier components:
    Assert(deconvIm);
    DVector* newRe = new DVector(C* (*deconvRe) - S*(*deconvIm));
    DVector* newIm = new DVector(S* (*deconvRe) + C*(*deconvIm));
    delete deconvRe;
    delete deconvIm;
    deconvRe = newRe;
    deconvIm = newIm;
  }

  // Make new derivatives if there are old ones
  if (ExposureGroup<T>::dRedxy) ExposureGroup<T>::createdXdY();

  // The variance vector for Fourier coeffs is unchanged by this op.
}

///////////////////////////////////////////////////////////////
// Derived class that does FT via native GL decomposition
///////////////////////////////////////////////////////////////

template <class T>
ExposureGroupGL<T>::ExposureGroupGL(PSFInformation& psf_,
				    int order,
				    double maskSigma_):
  b(order), maskSigma(maskSigma_),
  bReV(0), bImV(0),
  ExposureGroup<T>(psf_) {}

template <class T>
ExposureGroupGL<T>::~ExposureGroupGL() {
  if (bReV) delete bReV;
  if (bImV) delete bImV;
}


template <class T>
double
ExposureGroupGL<T>::prepare(Ellipse startBasis) {
  // Measure GL expansion of this galaxy
  const int initialOrder=4;	// Order of GL fit used to set basis
  GLSimple<T> gl(ExposureGroup<T>::felist, startBasis, initialOrder);
  gl.setCentering(false);	// Hold to common center of all images
  ExposureGroup<T>::phaseCenter = startBasis.getX0(); // remember FT origin
  gl.setMaskSigma(maskSigma);
  if (!gl.solve()) {
    // GL failure invalidates this image
    ExposureGroup<T>::flags = gl.getFlags();
    ExposureGroup<T>::flags |= GLFailure;
    return 0.;
  }

  Ellipse nativeBasis = gl.getBasis();

  // Get significance info:
  double f, varf;
  gl.b00(f, varf);

  // Conduct high-order GL fit at this basis:
  gl.setOrder(b.getOrder());
  b = gl.linearSolution();
  // Check for problems:
  ExposureGroup<T>::flags = gl.getFlags();
  if (gl.getFlags() & (Singularity + OutOfBounds + TooLarge + TooElliptical)) {
    ExposureGroup<T>::flags |= GLFailure;
    return 0;
  }

  // Make sure that target-domain k vectors and GL wavefunctions are set up
  static const int Nkpix = static_cast<int> (floor(kmaxGL/dkGL+0.01));


  // Build knorm vectors - first get all k points in the basis region
  vector<double> kxn, kyn;
  // The test points will cover only the right half-plane since
  // all functions are Hermitian
  // First half of the kx=0 line:
  for (int j=0; j<=Nkpix; j++) {
    kxn.push_back(0.);
    kyn.push_back(j*dkGL);
  }
  // Then all kx>0:
  for (int i=1; i<=Nkpix; i++) 
    for (int j=-Nkpix; j<=Nkpix; j++) {
      kxn.push_back(i*dkGL);
      kyn.push_back(j*dkGL);
    }

  // Get ktrue vectors and inverse PSF's for this native basis
  const int nk = kxn.size();
  if (ktruex) delete ktruex;
  if (ktruey) delete ktruey;
  if (deconvRe) delete deconvRe;
  if (deconvIm) delete deconvIm;

  // Make temporary vectors to hold useful points
  vector<double> tkx; tkx.reserve(nk);
  vector<double> tky; tky.reserve(nk);
  vector<double> tkxn; tkxn.reserve(nk);
  vector<double> tkyn; tkyn.reserve(nk);
  vector<double> invmtfRe; invmtfRe.reserve(nk);
  vector<double> invmtfIm; invmtfIm.reserve(nk);

  // Make a new Ellipse with same shear & mag as nativeE but no
  // displacement, so it can be used on k vectors:
  Ellipse nE = nativeBasis;
  nE.setX0(Position<double>(0.,0.));

  for (int i=0; i<nk; i++) {
    Position<double> knorm(kxn[i], kyn[i]);
    Position<double> ktrue = nE.inv(knorm);
    double kwx = ktrue.x; 
    double kwy = ktrue.y;

    DComplex mtf = ExposureGroup<T>::psf.T.kValue(Position<double> (kwx,kwy));
    // Place a floor on mtf to avoid division by zero.  Throw away
    // points outside this region and assume that T0sq will always
    // kill their influence.
    if (norm(mtf) < MTF_FLOOR*MTF_FLOOR) continue;

    tkx.push_back(kwx);
    tky.push_back(kwy);
    tkxn.push_back(kxn[i]);
    tkyn.push_back(kyn[i]);
    DComplex invmtf = 1./mtf;
    invmtfRe.push_back(invmtf.real());
    invmtfIm.push_back(invmtf.imag());
    // Mark origin as self-conjugate
    if (kwx==0. && kwy==0.) dcIndex = tkx.size()-1;
  }

  // Now have our culled set of points in k-space.
  ktruex = new DVector(tkx.size());
  std::copy(tkx.begin(), tkx.end(), ktruex->begin());
  ktruey = new DVector(tky.size());
  std::copy(tky.begin(), tky.end(), ktruey->begin());
  DVector knormx(tkxn.size());
  std::copy(tkxn.begin(), tkxn.end(), knormx.begin());
  DVector knormy(tkyn.size());
  std::copy(tkyn.begin(), tkyn.end(), knormy.begin());
  deconvRe = new DVector(tkx.size());
  deconvIm = new DVector(tkx.size());

  // Build new basis-function matrix - let LVector allocate them for us
  DMatrix* BRe = 0;
  DMatrix* BIm = 0;
  LVector::kBasis(knormx, knormy,
		  BRe, BIm,
		  b.getOrder());
  // Multiply the DC frequency point by 0.5 since we are
  // always going to integrate over half-planes
  BRe->row(dcIndex) *= 0.5;
  // divide basis matrices by MTF
  // Add scaling to inverse-MTF so results of dot products will look like integrals:
  tmv::DiagMatrix<double> invMTFRe(invmtfRe.size());
  std::copy(invmtfRe.begin(), invmtfRe.end(), invMTFRe.begin());
  tmv::DiagMatrix<double> invMTFIm(invmtfIm.size());
  std::copy(invmtfIm.begin(), invmtfIm.end(), invMTFIm.begin());

  // ??? scaling of these ???
  double nEdet = nE.getMatrix().det();
  ExposureGroup<T>::ftnorm = 2. * dkGL * dkGL / nEdet;

  invMTFRe *= 2*dkGL*dkGL / nEdet;  // ??? check on this last term ???
  invMTFIm *= 2*dkGL*dkGL / nEdet;

  DMatrix basisRe =  invMTFRe * (*BRe) - invMTFIm * (*BIm);
  DMatrix basisIm =  invMTFRe * (*BIm) + invMTFIm * (*BRe);
  delete BRe; BRe = 0;
  delete BIm; BIm = 0;

  // Project b vector into Fourier coefficients:
  Assert(basisRe.ncols() == b.rVector().size());
  *deconvRe = basisRe * b.rVector();
  *deconvIm = basisIm * b.rVector();

  // save away appropriate covariance matrices.
  // If the inverse-covariance of the b coefficients is already
  // known from GL fit, and has an SVD C = V^T S V, then it is most efficient
  // to save away the matrices
  //  basis * V^T * (1/sqrt(S)) 
  DVector invsqrtS(gl.fisherSV_S());
  double maxS;
  double minS;
  for (int i=0; i<invsqrtS.size(); i++) {
    if (i==0) {maxS=invsqrtS[i]; minS=invsqrtS[i];}
    maxS = MAX(maxS, invsqrtS[i]);
    minS = MIN(minS, invsqrtS[i]);
    invsqrtS[i] = 1./sqrt(invsqrtS[i]);
  }
  // Set flag for (near) singularities here
  if (maxS/minS > ILL_CONDITION) ExposureGroup<T>::flags |= BDegeneracy;

  DMatrix tmp2 = gl.fisherSV_V().transpose() * tmv::DiagMatrix<double>(invsqrtS);
  if (bReV) delete bReV;
  bReV = new DMatrix(basisRe * tmp2);
  if (bImV) delete bImV;
  bImV = new DMatrix(basisIm * tmp2);

  ExposureGroup<T>::createdXdY();

  return f/sqrt(varf);
}

template <class T>
tmv::SymMatrix<double>
ExposureGroupGL<T>::covarianceOf(const DMatrix& mRe, 
				 const DMatrix& mIm) {
  throw FDNTError("ExposureGroupGL<T>::covarianceOf() not yet implemented");

  const int nRe = mRe.nrows();
  const int nIm = mIm.nrows();
  const int nTests = nRe + nIm;
  tmv::SymMatrix<double> cov(nTests);

  DMatrix leftRe = mRe * (*bReV);
  DMatrix leftIm = mIm * (*bReV);

  // And we can construct the covariance matrix of these:
  // Use TMV special function when we know output is symmetric:
  if (nRe>0) {
    tmv::SymMultMM<false>(1.,
			  leftRe.view(),
			  leftRe.transpose(),
			  cov.subSymMatrix(0,nRe));
  }
  if (nIm > 0) {
    tmv::SymMultMM<false>(1.,
			  leftIm.view(),
			  leftIm.transpose(),
			  cov.subSymMatrix(nRe,nTests));
    if (nRe>0) {
      // Cross terms:
      cov.subMatrix(0,nRe,nRe,nTests) = leftRe * leftIm.transpose();
    }
  }
  return cov;
}

template <class T>
void
ExposureGroupGL<T>::setPhaseCenter(double x0, double y0) {
  double dx = x0 - ExposureGroup<T>::phaseCenter.x;
  double dy = y0 - ExposureGroup<T>::phaseCenter.y;
  ExposureGroup<T>::phaseCenter.x = x0;
  ExposureGroup<T>::phaseCenter.y = y0;
  if (! (deconvRe || deconvIm || ExposureGroup<T>::dRedxy || ExposureGroup<T>::dImdxy 
	 || bReV || bImV) )
    return;	// Nothing to change!
  // Build cos & sin vectors as diag matrices
  Assert(ktruex);
  const int nk=ktruex->size();
  tmv::DiagMatrix<double> C(nk, 0.);
  tmv::DiagMatrix<double> S(nk, 0.);
  DVector phase = dx * (*ktruex) + dy * (*ktruey);
  for (int i=0; i<nk; i++) {
    C(i) = cos(phase[i]);
    S(i) = sin(phase[i]);
  }
  if (deconvRe) {
    // Phase shift the Fourier components:
    Assert(deconvIm);
    DVector* newRe = new DVector(C* (*deconvRe) - S*(*deconvIm));
    DVector* newIm = new DVector(S* (*deconvRe) + C*(*deconvIm));
    delete deconvRe;
    delete deconvIm;
    deconvRe = newRe;
    deconvIm = newIm;
  }

  // Make new derivatives if there are old ones
  if (ExposureGroup<T>::dRedxy) ExposureGroup<T>::createdXdY();

  if (bReV) {
    // Create new covariance matrices if there are old ones
    Assert(bImV);
    DMatrix* newbReV = new DMatrix( C * (*bReV) - S*(*bImV));
    DMatrix* newbImV = new DMatrix( S * (*bReV) + C*(*bImV));
    delete bReV;
    delete bImV;
    bReV = newbReV;
    bImV = newbImV;
  }
}

/////////////////////
// Instantiate:
/////////////////////
template class ExposureGroup<float>;
template class ExposureGroupGL<float>;
template class ExposureGroupFT<float>;
