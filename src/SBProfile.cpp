// $Id: SBProfile.cpp,v 1.17 2012/01/11 14:18:51 dgru Exp $
//
// Functions for the Surface Brightness Profile Class
//
#include "SBProfile.h"
#include <cmath>
#include "Simpson.h"
#include "UseTMV.h"
#include "Solve.h"

using namespace sbp;

// ????? Change treatement of aliased images to simply add in the aliased
// FT components instead of doing a larger FT and then subsampling!
// ??? Make a formula for asymptotic high-k SBSersic::kValue ??

// Parameters controlling behavior of all classes:
const int 
SBProfile::MINIMUM_FFT_SIZE = 128;
const int 
SBProfile::MAXIMUM_FFT_SIZE = 4096;
// Allow aliasing of Fourier modes below this amplitude, roughly.
// Also set the FFT image size such that this fraction of flux (or less) is "folded."
const double 
SBProfile::ALIAS_THRESHOLD = 0.005;	


//
// Virtual methods of Base Class "SBProfile"
//

SBProfile* 
SBProfile::distort(const Ellipse e) const {
  return new SBDistort(*this,e);
}

SBProfile* 
SBProfile::rotate(const double theta) const {
  return new SBDistort(*this,cos(theta),-sin(theta),sin(theta),cos(theta));
}

SBProfile* 
SBProfile::shift(double dx, double dy) const {
  return new SBDistort(*this,1.,0.,0.,1., Position<double>(dx,dy));
}

//
// Common methods of Base Class "SBProfile"
//

Image<> 
SBProfile::draw(double dx, int wmult) const {
  Image<> img;
  draw(img, dx, wmult);
  return img;
}

double
SBProfile::draw(Image<> img, double dx, int wmult) const {
  if (isAnalyticX())
    return plainDraw(img, dx, wmult);
  else
    return fourierDraw(img, dx, wmult);
}

// First is a simple case wherein we have a formula for x values:
double
SBProfile::plainDraw(Image<> I, double dx, int wmult) const {
  // Determine desired dx:
  if (dx<=0.) dx = PI / maxK();
  if (!I.getBounds()) {
    if (wmult<1) throw SBError("Requested wmult<1 in plainDraw()");
    // Need to choose an image size
    int N = static_cast<int> (ceil(2*PI/(dx*stepK())));

    // Round up to an even value
    N = 2*( (N+1)/2);
    N *= wmult;	// make even bigger if desired
    Bounds<int> imgsize(-N/2, N/2-1, -N/2, N/2-1);
    I.resize(imgsize);
  } else {
    // recenter an existing image, to be consistent with fourierDraw:
    int xSize = I.XMax()-I.XMin()+1, ySize = I.YMax()-I.YMin()+1;
    I.shift(-xSize/2, -ySize/2);
  }

  return fillXImage(I, dx);
}

double
SBProfile::fillXImage(Image<> I, double dx) const {
  double totalflux=0;
  for (int y = I.YMin(); y <= I.YMax(); y++) {
    int x = I.XMin(); 
    Image<>::iter ee=I.rowEnd(y);
    for (Image<>::iter it=I.rowBegin(y);
	 it!=ee;
	 ++it, ++x) {
      Position<double> p(x*dx,y*dx); // since x,y are pixel indices
      *it = xValue(p);
#ifdef DANIELS_TRACING
      cout << "x=" << x << ", y=" << y << ": " << *it << endl;
      cout << "--------------------------" << endl;
#endif
      totalflux += *it;
    } 
  }
  // add image header to include "dx" info
  try {I.header()->erase("DK");} catch (ImageHeaderError &i) {}
  I.header()->replace("DX", dx);
  return totalflux * (dx*dx);
}

// Now the more complex case: real space via FT from k space.
// Will enforce image size is power of 2 or 3x2^n.
// Aliasing will be handled by folding the k values before transforming
// And enforce no image folding
/**/
/**/
//**/ #define DEBUG
double
SBProfile::fourierDraw(Image<> I, double dx, int wmult) const {
  Bounds<int> imgBounds;	// Bounds for output image
  bool sizeIsFree = !I.getBounds();
  if (wmult<1) throw SBError("Requested wmult<1 in fourierDraw()");
  // First choose desired dx if we were not given one:
  if (dx<=0.) {
    // Choose for ourselves:
    dx = PI / maxK();
  }

#ifdef DEBUG
  cerr << " maxK() " << maxK() << " dx " << dx << endl;
#endif

  // Now decide how big the FT must be to avoid folding:
  double xRange = 2*PI*wmult / stepK();
  int Nnofold = static_cast<int> (ceil(xRange / dx -0.0001));
#ifdef DEBUG
  cerr << " stepK() " << stepK() << " Nnofold " << Nnofold << endl;
#endif
  
  // And if there is a target image size, we must make something big enough to cover
  // the target image size:
  if (!sizeIsFree) {
    int xSize, ySize;
    xSize = I.XMax()-I.XMin()+1;
    ySize = I.YMax()-I.YMin()+1;
    if (xSize  > Nnofold) Nnofold = xSize;
    if (ySize  > Nnofold) Nnofold = ySize;
    xRange = Nnofold * dx;
  }

  // Round up to a good size for making FFTs:
  int NFT = fft::goodFFTSize(Nnofold);
  NFT = MAX(NFT,MINIMUM_FFT_SIZE);
#ifdef DEBUG
  cerr << " After adjustments: Nnofold " << Nnofold << " NFT " << NFT << endl;
#endif
  if (NFT > MAXIMUM_FFT_SIZE)
    FormatAndThrow<SBError>() << "fourierDraw() requires an FFT that is too large, " << NFT;

  // If we are free to set up output image, make it size of FFT
  if (sizeIsFree) {
    int Nimg = NFT;
    // Reduce to make even
    Nimg = 2*(Nimg/2);
    imgBounds = Bounds<int>(-Nimg/2, Nimg/2-1, -Nimg/2, Nimg/2-1);
    I.resize(imgBounds);
  } else {
    // Going to shift the output image to be centered near zero
    int xSize, ySize;
    xSize = I.XMax()-I.XMin()+1;
    ySize = I.YMax()-I.YMin()+1;
    I.shift(-xSize/2, -ySize/2);
  }
  double dk = 2.*PI/(NFT*dx);
#ifdef DEBUG
  cerr << " After adjustments: dx " << dx << " dk " << dk << " maxK " << dk*NFT/2 << endl;
#endif
  Assert(dk <= stepK());
  XTable* xtmp=0;
  if (NFT*dk/2 > maxK()) {
    // No aliasing: build KTable and transform
    KTable kt(NFT,dk);
    //cerr << "no aliasing, NFT " << NFT << " dk " << dk << endl;
    fillKGrid(kt); 
    xtmp = kt.transform();
  } else {
    // There will be aliasing.  Construct a KTable out to maxK() and
    // then wrap it
    int Nk = static_cast<int> (ceil(maxK()/dk)) * 2;
    KTable kt(Nk, dk);
    fillKGrid(kt);
    //cerr << "aliasing, NFT " << NFT << " dk " << dk << " Nk " << Nk << endl;
    KTable* kt2 = kt.wrap(NFT);
    xtmp = kt2->transform();
    delete kt2;
  }
  int Nxt = xtmp->getN();
  Bounds<int> xb(-Nxt/2, Nxt/2-1, -Nxt/2, Nxt/2-1);
  if (I.YMin() < xb.getYMin()
      || I.YMax() > xb.getYMax()
      || I.XMin() < xb.getXMin()
      || I.XMax() > xb.getXMax()) {
    /**/cerr << "Bounds error!! target image bounds " << I.getBounds()
	     << " and FFT range " << xb
	     << endl;
    throw SBError("fourierDraw() FT bounds do not cover target image");
  }
  double sum=0.;
  for (int y = I.YMin(); y <= I.YMax(); y++)
    for (int x = I.XMin(); x <= I.XMax(); x++) {
      I(x,y) = xtmp->xval(x,y);
      sum += I(x,y);
    }

  try {I.header()->erase("DK");} catch (ImageHeaderError &i) {}
  I.header()->replace("DX", dx);

  delete xtmp;  // no memory leak!
  return sum*dx*dx;;
}


void
SBProfile::drawK(Image<> Re, Image<> Im, double dk, int wmult) const {
  if (isAnalyticK()) 
    plainDrawK(Re, Im, dk, wmult);   // calculate in k space
  else               
    fourierDrawK(Re, Im, dk, wmult); // calculate via FT from real space
  return;
}

void
SBProfile::plainDrawK(Image<> Re, Image<> Im, double dk, int wmult) const {
  // Make sure input images match or are both null
  Assert(!(Re.getBounds() || Im.getBounds()) || (Re.getBounds() == Im.getBounds()));
  if (dk<=0.) dk = stepK();
  
  if (!Re.getBounds()) {
    if (wmult<1) throw SBError("Requested wmult<1 in plainDrawK()");
    // Need to choose an image size
    int N = static_cast<int> (ceil(2.*maxK()*wmult / dk));
    // Round up to an even value
    N = 2*( (N+1)/2);

    Bounds<int> imgsize(-N/2, N/2-1, -N/2, N/2-1);
    Re.resize(imgsize);
    Im.resize(imgsize);
  } else {
    // recenter an existing image, to be consistent with fourierDrawK:
    int xSize = Re.XMax()-Re.XMin()+1, ySize = Re.YMax()-Re.YMin()+1;
    Re.shift(-xSize/2, -ySize/2);
    Im.shift(-xSize/2, -ySize/2);
  }

  // ??? Make this into a virtual function to allow pipelining?
  for (int y = Re.YMin(); y <= Re.YMax(); y++) {
    int x = Re.XMin(); 
    Image<>::iter ee=Re.rowEnd(y);
    Image<>::iter it;
    Image<>::iter it2;
    for (it=Re.rowBegin(y), it2=Im.rowBegin(y);
	 it!=ee;
	 ++it, ++it2, ++x) {
      Position<double> p(x*dk,y*dk); // since x,y are pixel indicies
      DComplex c = this->kValue(p);  
      *it = c.real(); 
      *it2 = c.imag(); 
    } 
  }

  // add image header to include dk
  try {Re.header()->erase("DX");} catch (ImageHeaderError &i) {}
  Re.header()->replace("DK", dk);

  try {Im.header()->erase("DX");} catch (ImageHeaderError &i) {}
  Im.header()->replace("DK", dk);

  return;
}

// Build K domain by transform from X domain.  This is likely
// to be a rare event but what the heck.  Enforce no "aliasing"
// by oversampling and extending x domain if needed.  Force
// power of 2 for transform

void
SBProfile::fourierDrawK(Image<> Re, Image<> Im, double dk, int wmult) const {
  Assert(!(Re.getBounds() || Im.getBounds()) || (Re.getBounds() == Im.getBounds()));

  int oversamp =1;	// oversampling factor
  Bounds<int> imgBounds;	// Bounds for output image
  bool sizeIsFree = !Re.getBounds();
  if (wmult<1) throw SBError("Requested wmult<1 in fourierDrawK()");
  bool canReduceDk=true;
  // First choose desired dx
  if (dk<=0.) {
    // Choose for ourselves:
    dk = stepK();
    canReduceDk = true;
  } else {
    // We have a value we must produce.  Do we need to oversample in k
    // to avoid folding from real space?
    // Note a little room for numerical slop before triggering oversampling:
    oversamp = static_cast<int> ( ceil(dk/stepK() - 0.0001));
    canReduceDk = false;	// Force output image to input dx.
  }

  // Now decide how big the FT must be to avoid folding
  double kRange = 2*maxK()*wmult;
  int Nnofold = static_cast<int> (ceil(oversamp*kRange / dk -0.0001));
  
  // And if there is a target image size, we must make something big enough to cover
  // the target image size:
  if (!sizeIsFree) {
    int xSize, ySize;
    xSize = Re.XMax()-Re.XMin()+1;
    ySize = Re.YMax()-Re.YMin()+1;
    if (xSize * oversamp > Nnofold) Nnofold = xSize*oversamp;
    if (ySize * oversamp > Nnofold) Nnofold = ySize*oversamp;
    kRange = Nnofold * dk / oversamp;
    // If the input image *size* was specified but not the input *dk*, then
    // we will hold dk at the Nyquist scale:
    canReduceDk = false;
  }

  // Round up to a power of 2 to get required FFT size
  int NFT = MINIMUM_FFT_SIZE;
  while (NFT < Nnofold && NFT<= MAXIMUM_FFT_SIZE) NFT *= 2;
  if (NFT > MAXIMUM_FFT_SIZE)
    throw SBError("fourierDrawK() requires an FFT that is too large");

  // If we are free to set up output image, make it size of FFT less oversampling
  if (sizeIsFree) {
    int Nimg = NFT / oversamp;
    // Reduce to make even
    Nimg = 2*(Nimg/2);
    imgBounds = Bounds<int>(-Nimg/2, Nimg/2-1, -Nimg/2, Nimg/2-1);
    Re.resize(imgBounds);
    Im.resize(imgBounds);
    // Reduce dk if 2^N made left room to do so.
    if (canReduceDk) {
      dk = kRange / Nimg; 
    }
  } else {
    // Going to shift the output image to be centered near zero
    int xSize, ySize;
    xSize = Re.XMax()-Re.XMin()+1;
    ySize = Re.YMax()-Re.YMin()+1;
    Re.shift(-xSize/2, -ySize/2);
    Im.shift(-xSize/2, -ySize/2);
  }

  double dx = 2.*PI*oversamp/(NFT*dk);
  XTable xt(NFT,dx);
  this->fillXGrid(xt);
  KTable *ktmp = xt.transform();

  int Nkt = ktmp->getN();
  Bounds<int> kb(-Nkt/2, Nkt/2-1, -Nkt/2, Nkt/2-1);
  if (Re.YMin() < kb.getYMin()
      || Re.YMax()*oversamp > kb.getYMax()
      || Re.XMin()*oversamp < kb.getXMin()
      || Re.XMax()*oversamp > kb.getXMax()) {
    /**/cerr << "Bounds error!! oversamp is " << oversamp
	     << " target image bounds " << Re.getBounds()
	     << " and FFT range " << kb
	     << endl;
    throw SBError("fourierDrawK() FT bounds do not cover target image");
  }

  for (int y = Re.YMin(); y <= Re.YMax(); y++)
    for (int x = Re.XMin(); x <= Re.XMax(); x++) {
      Re(x,y) = ktmp->kval(x*oversamp,y*oversamp).real();
      Im(x,y) = ktmp->kval(x*oversamp,y*oversamp).imag();
    }

  try {Re.header()->erase("DX");} catch (ImageHeaderError &i) {}
  Re.header()->replace("DK", dk);
  try {Im.header()->erase("DX");} catch (ImageHeaderError &i) {}
  Im.header()->replace("DK", dk);

  delete ktmp;  // no memory leak!
}

void
SBProfile::fillXGrid(XTable &xt) const {
  int N = xt.getN();
  double dx = xt.getDx();
  for (int iy = -N/2; iy < N/2; iy++)
    for (int ix = -N/2; ix < N/2; ix++) {
      Position<double> x(ix*dx,iy*dx);
      xt.xSet(ix,iy,xValue(x));
    }
  return;
}

void
SBProfile::fillKGrid(KTable &kt) const {
  int N = kt.getN();
  double dk = kt.getDk();
  for (int iy = -N/2; iy < N/2; iy++) {
    // Only need ix>=0 because it's Hermitian:
    for (int ix = 0; ix <= N/2; ix++) {
      Position<double> k(ix*dk,iy*dk);
      kt.kSet(ix,iy,kValue(k));
    }
  }
  return;
}

//
// Methods for Derived Classes
//

void
SBAdd::initialize() {
  sumflux = sumfx = sumfy = 0.;
  maxMaxK = minStepK = 0.;
  allAxisymmetric = allAnalyticX = allAnalyticK = true;
}

void
SBAdd::add(const SBProfile& rhs, double scale) {
  // Need a non-const copy of the rhs:
  SBProfile* p=rhs.duplicate();
  

  // Keep track of where first new summand is on list:
  list<SBProfile*>::iterator newptr = plist.end();

  // Add new summand(s) to the plist:
  SBAdd *sba = dynamic_cast<SBAdd*> (p);
  if (sba) {  
    // If rhs is an SBAdd, copy its full list here
    list<SBProfile*>::const_iterator pptr;
    for (pptr = sba->plist.begin(); pptr!=sba->plist.end(); ++pptr) {
      if (newptr==plist.end()) {
	plist.push_back((*pptr)->duplicate()); 
	// Rescale flux for duplicate copy if desired:
	if (scale!=1.) 
	  plist.back()->setFlux( scale*plist.back()->getFlux());
	newptr = --plist.end();  // That was first new summand
      } else {
	plist.push_back((*pptr)->duplicate()); 
      }
    }
    delete sba; // no memory leak! 
  } else {
    plist.push_back(p);
    // Rescale flux for duplicate copy if desired:
    if (scale!=1.) 
      plist.back()->setFlux( scale*plist.back()->getFlux());
    newptr = --plist.end();  // That was first new summand
  }

  // Accumulate properties of all summands
  while (newptr != plist.end()) {
    sumflux += (*newptr)->getFlux();
    sumfx += (*newptr)->getFlux() * (*newptr)->centroidX();
    sumfy += (*newptr)->getFlux() * (*newptr)->centroidY();
    if ( (*newptr)->maxK() > maxMaxK) maxMaxK = (*newptr)->maxK();
    if ( minStepK<=0. || ((*newptr)->stepK() < minStepK)) minStepK = (*newptr)->stepK();
    allAxisymmetric = allAxisymmetric && (*newptr)->isAxisymmetric();
    allAnalyticX = allAnalyticX && (*newptr)->isAnalyticX();
    allAnalyticK = allAnalyticK && (*newptr)->isAnalyticK();
    newptr++;
  }
  return; 
}

double 
SBAdd::xValue(Position<double> _p) const {
  double xv = 0.;  
  list<SBProfile*>::const_iterator pptr;
  for (pptr = plist.begin(); pptr != plist.end(); ++pptr)
  {
    xv += (*pptr)->xValue(_p);
  }
  return xv;
} 

DComplex 
SBAdd::kValue(Position<double> _p) const {
  DComplex kv = 0.;  
  list<SBProfile*>::const_iterator pptr;
  for (pptr = plist.begin(); pptr != plist.end(); ++pptr)
    kv += (*pptr)->kValue(_p);
  return kv;
} 

void
SBAdd::fillKGrid(KTable &kt) const {
  if (plist.empty()) kt.clear();
  list<SBProfile*>::const_iterator pptr = plist.begin();
  (*pptr)->fillKGrid(kt);
  ++pptr;
  if (pptr==plist.end()) return;
  int N = kt.getN();
  double dk = kt.getDk();
  KTable k2(N,dk);
  for ( ; pptr!= plist.end(); ++pptr) {
    (*pptr)->fillKGrid(k2);
    kt.accumulate(k2);
  }
}

void
SBAdd::fillXGrid(XTable &xt) const {
  if (plist.empty()) xt.clear();
  list<SBProfile*>::const_iterator pptr = plist.begin();
  (*pptr)->fillXGrid(xt);
  ++pptr;
  if (pptr==plist.end()) return;
  int N = xt.getN();
  double dx = xt.getDx();
  XTable x2(N,dx);
  for ( ; pptr!= plist.end(); ++pptr) {
    (*pptr)->fillXGrid(x2);
    xt.accumulate(x2);
  }
}

void
SBAdd::setCentroid(Position<double> p_) {
  // Try passing this on to all summands.  If we were axisymmetric before
  // we won't be now...
  allAxisymmetric = false;
  Position<double> d1 = p_ - centroid(), d2;
  sumfx += sumflux*d1.x;
  sumfy += sumflux*d1.y;
  list<SBProfile*>::const_iterator pptr; 
  for (pptr = plist.begin(); pptr != plist.end(); ++pptr) {
    d2 = (*pptr)->centroid() + d1;
    (*pptr)->setCentroid(d2);
  }
  return;
}

void
SBAdd::setFlux(double f) {
  if (sumflux==0.) throw SBError("SBAdd::setFlux not possible when flux=0 to start");
  double m = f/getFlux();  // Factor by which to change flux
  list<SBProfile*>::iterator pptr; 
  for (pptr = plist.begin(); pptr != plist.end(); ++pptr) {
    double pf = (*pptr)->getFlux();  
    (*pptr)->setFlux(pf*m);
  }
  sumflux *=m;
  sumfx *= m;
  sumfy *= m;
  return;
}

//
// "SBDistort" Class 
//
SBDistort::SBDistort(const SBProfile& sbin, 
		     double mA, double mB, double mC, double mD,
		     Position<double> x0_):
  matrixA(mA), matrixB(mB), matrixC(mC), matrixD(mD), x0(x0_)
{
  SBProfile* p=sbin.duplicate();
  SBDistort* sbd = dynamic_cast<SBDistort*> (p);
  if (sbd) {
    // We are distorting something that's already a distortion.
    // So just compound the affine transformaions
    adaptee = sbd->adaptee->duplicate();
    x0 = x0_ + fwd(sbd->x0);
    // New matrix is product (M_this) * (M_old)
    matrixA = mA*sbd->matrixA + mB*sbd->matrixC;
    matrixB = mA*sbd->matrixB + mB*sbd->matrixD;
    matrixC = mC*sbd->matrixA + mD*sbd->matrixC;
    matrixD = mC*sbd->matrixB + mD*sbd->matrixD;
    delete sbd;
  } else {
    // Distorting something generic
    adaptee = p;
  }
  initialize();
}
    
SBDistort::SBDistort(const SBProfile& sbin, const Ellipse e_) {
  // First get what we need from the Ellipse:
  DMatrix m = e_.getMatrix();
  matrixA = m(0,0);
  matrixB = m(0,1);
  matrixC = m(1,0);
  matrixD = m(1,1);
  x0 = e_.getX0();
  // Then repeat generic construction:
  SBProfile* p=sbin.duplicate();
  SBDistort* sbd = dynamic_cast<SBDistort*> (p);
  if (sbd) {
    // We are distorting something that's already a distortion.
    // So just compound the affine transformaions
    adaptee = sbd->adaptee->duplicate();
    x0 = e_.getX0() + fwd(sbd->x0);
    // New matrix is product (M_this) * (M_old)
    double mA = matrixA; double mB=matrixB; double mC=matrixC; double mD=matrixD;
    matrixA = mA*sbd->matrixA + mB*sbd->matrixC;
    matrixB = mA*sbd->matrixB + mB*sbd->matrixD;
    matrixC = mC*sbd->matrixA + mD*sbd->matrixC;
    matrixD = mC*sbd->matrixB + mD*sbd->matrixD;
    delete sbd;
  } else {
    // Distorting something generic
    adaptee = p;
  }
  initialize();
}

void
SBDistort::initialize() {
  double det = matrixA*matrixD-matrixB*matrixC;
  if (det==0.) throw SBError("Attempt to SBDistort with degenerate matrix");
  absdet = abs(det);
  invdet = 1./det;

  //**/cerr << "Matrix: " << matrixA << " " << matrixB
  //	   << " " << matrixC << " " << matrixD << endl;
  double h1 = hypot( matrixA+matrixD, matrixB-matrixC);
  double h2 = hypot( matrixA-matrixD, matrixB+matrixC);
  major = 0.5*abs(h1+h2);
  minor = 0.5*abs(h1-h2);
  if (major<minor) SWAP(major,minor);
  stillIsAxisymmetric = adaptee->isAxisymmetric() 
    && (matrixB==-matrixC) 
    && (matrixA==matrixD)
    && (x0.x==0.) && (x0.y==0.);	// Need pure rotation
  //**/cerr << "Set major = " << major << " and minor " << minor << endl;
}

// Specialization of fillKGrid is desired since the phase terms from shift 
// are factorizable:
void
SBDistort::fillKGrid(KTable &kt) const {
  double N = (double) kt.getN();
  double dk = kt.getDk();

  if (x0.x==0. && x0.y==0.) {
    // Branch to faster calculation if there is no centroid shift:
    for (int iy = -N/2; iy < N/2; iy++) {
      // only need ix>=0 since it's Hermitian:
      for (int ix = 0; ix <= N/2; ix++) {
	Position<double> k(ix*dk,iy*dk);
	kt.kSet(ix,iy,kValNoPhase(k));
      }
    }
  } else {
    DComplex dxexp(0,-dk*x0.x),   dyexp(0,-dk*x0.y);
    DComplex dxphase(exp(dxexp)), dyphase(exp(dyexp));
    // xphase, yphase: current phase value
    DComplex yphase(exp(-dyexp*N/2.));
    for (int iy = -N/2; iy < N/2; iy++) {
      DComplex phase = yphase; // since kx=0 to start
      // Only ix>=0 since it's Hermitian:
      for (int ix = 0; ix <= N/2; ix++) {
	Position<double> k(ix*dk,iy*dk);
	kt.kSet(ix,iy,kValNoPhase(k) * phase);
	phase *= dxphase;
      }
      yphase *= dyphase;
    }
  }
  return;
}

//
// SBConvolve class - adding new members
//
void 
SBConvolve::add(const SBProfile& rhs) {
  if (!rhs.isAnalyticK()) throw SBError("SBConvolve requires members to be analytic in k");
  // If this is the first thing being added to the list, initialize some accumulators
  if (plist.empty()) {
    x0 = y0 = 0.;
    fluxProduct = 1.;
    minMaxK = 0.;
    minStepK = 0.;
    isStillAxisymmetric = true;
  }

  // Need a non-const copy of the rhs:
  SBProfile* p=rhs.duplicate();

  // Keep track of where first new term is on list:
  list<SBProfile*>::iterator newptr = plist.end();

  // Add new terms(s) to the plist:
  SBConvolve *sbc = dynamic_cast<SBConvolve*> (p);
  if (sbc) {  
    // If rhs is an SBConvolve, copy its list here
    fluxScale *= sbc->fluxScale;
    list<SBProfile*>::iterator pptr;
    for (pptr = sbc->plist.begin(); pptr!=sbc->plist.end(); ++pptr) {
      if (newptr==plist.end()) {
	plist.push_back((*pptr)->duplicate()); 
	newptr = --plist.end();  // That was first new term
      } else {
	plist.push_back((*pptr)->duplicate()); 
      }
    }
    delete sbc; // no memory leak! 
  } else {
    plist.push_back(p);
    newptr = --plist.end();  // That was first new term
  }

  // Accumulate properties of all terms
  while (newptr != plist.end()) {
    fluxProduct *= (*newptr)->getFlux();
    x0 += (*newptr)->centroidX();
    y0 += (*newptr)->centroidY();
    if ( minMaxK<=0. || (*newptr)->maxK() < minMaxK) minMaxK = (*newptr)->maxK();
    if ( minStepK<=0. || ((*newptr)->stepK() < minStepK)) minStepK = (*newptr)->stepK();
    isStillAxisymmetric = isStillAxisymmetric && (*newptr)->isAxisymmetric();
    newptr++;
  }
}

void
SBConvolve::setCentroid(Position<double> p_) {
  // Try passing this on to the first of the convolvees.  Often will throw...
  // Also any previous axisymmetry will be gone:
  if (plist.empty()) return;	// nothing to do
  isStillAxisymmetric = false;
  Position<double> newFront = p_ - centroid() + plist.front()->centroid();
  plist.front()->setCentroid(newFront);
  x0 = p_.x; y0 = p_.y;
}

void
SBConvolve::fillKGrid(KTable &kt) const {
  if (plist.empty()) kt.clear();
  list<SBProfile*>::const_iterator pptr = plist.begin();
  (*pptr)->fillKGrid(kt);
  kt *= fluxScale;
  ++pptr;
  if (pptr==plist.end()) return;
  int N = kt.getN();
  double dk = kt.getDk();
  KTable k2(N,dk);
  for ( ; pptr!= plist.end(); ++pptr) {
    (*pptr)->fillKGrid(k2);
    kt*=k2;
  }
}

//
// "SBGaussian" Class 
//
double
SBGaussian::xValue(Position<double> p) const{
  double r2 = p.x*p.x + p.y*p.y;
  double xval = flux * exp( -r2/2./sigma/sigma );
  xval /= 2*PI*sigma*sigma;  // normalize
  return xval;
}

DComplex
SBGaussian::kValue(Position<double> p) const{
  double r2 = p.x*p.x + p.y*p.y;
  DComplex kval(flux * exp(-(r2)*sigma*sigma/2.),0);
  return kval;
}


//
// SBExponential Class
//

// Set stepK so that folding occurs when excluded flux=ALIAS_THRESHOLD
// Or at least 6 scale lengths
double
SBExponential::stepK() const{
  // A fast solution to (1+R)exp(-R)=ALIAS_THRESHOLD:
  double R = -log(ALIAS_THRESHOLD);
  for (int i=0; i<3; i++) R = -log(ALIAS_THRESHOLD) + log(1+R);
  R = MAX(6., R);
  return PI / (R*r0);
}

double
SBExponential::xValue(Position<double> p) const{
  double r = sqrt(p.x*p.x + p.y*p.y);
  double xval = flux * exp(-r/r0);
  xval /= r0*r0*2*PI;   // normalize
  return xval;
}

DComplex
SBExponential::kValue(Position<double> p) const {
  double kk = p.x*p.x+p.y*p.y;
  double temp = 1 + kk*r0*r0;         // [1+k^2*r0^2]
  DComplex kval( flux/sqrt(temp*temp*temp), 0.);
  return kval;
}

//
// SBAiry Class
//

// Note x & y are in units of lambda/D here.  Integral over area
// will give unity in this normalization.

double
SBAiry::xValue(Position<double> p) const {
  double radius = sqrt(p.x*p.x+p.y*p.y);
  double nu = radius*PI*D;
  double xval;
  if (nu<0.01)
    // lim j1(u)/u = 1/2
    xval =  D * (1-obscuration*obscuration);
  else {
    xval = 2*D*( j1(nu) - obscuration*j1(obscuration*nu)) /
      nu ;	//See Schroeder eq (10.1.10)
  }
  xval*=xval;
  // Normalize to give unit flux integrated over area.
  xval /= (1-obscuration*obscuration)*4./PI;
  return xval*flux;
}

double 
SBAiry::chord(const double r, const double h) const {
  if (r<h) throw SBError("Airy calculation r<h");
  else if (r==0.) return 0.;
  else if (r<0 || h<0) throw SBError("Airy calculation (r||h)<0");
  return r*r*asin(h/r) -h*sqrt(r*r-h*h);
}

/* area inside intersection of 2 circles radii r & s, seperated by t*/
double 
SBAiry::circle_intersection(double r, double s, double t) const {
  double h;
  if (r<0. || s<0.) throw SBError("Airy calculation negative radius");
  t = fabs(t);
  if (t>= r+s) return 0.;
  if (r<s) {
    double temp;
    temp = s;
    s = r;
    r = temp;
  }
  if (t<= r-s) return PI*s*s;

  /* in between we calculate half-height at intersection */
  h = 0.5*(r*r + s*s) - (pow(t,4.) + (r*r-s*s)*(r*r-s*s))/(4.*t*t);
  if (h<0) {
    throw SBError("Airy calculation half-height invalid");
  }
  h = sqrt(h);

  if (t*t < r*r - s*s) 
    return PI*s*s - chord(s,h) + chord(r,h);
  else
    return chord(s,h) + chord(r,h);
}

/* area of two intersecting identical annuli */
double
SBAiry::annuli_intersect(double r1, double r2, double t) const {
  if (r1<r2) {
    double temp;
    temp = r2;
    r2 = r1;
    r1 = temp;
  }
  return circle_intersection(r1,r1,t)
    - 2 * circle_intersection(r1,r2,t)
    +  circle_intersection(r2,r2,t);
}

/* Beam pattern of annular aperture, in k space, which is just the
 * autocorrelation of two annuli.  Normalize to unity at k=0 for now */
double
SBAiry::annuli_autocorrelation(const double k) const {
  double k_scaled = k / (PI*D);
  double norm = PI*(1. - obscuration*obscuration);
  return annuli_intersect(1.,obscuration,k_scaled)/norm;
}

DComplex
SBAiry::kValue(Position<double> p) const{
  double radius = sqrt(p.x*p.x+p.y*p.y);
  // calculate circular FT(PSF) on p'=(x',y')
  double r = annuli_autocorrelation(radius);
  DComplex kval(r, 0.);
  return kval*flux;
}


//
// SBBox Class
//

double
SBBox::xValue(Position<double> p) const {
  if (fabs(p.x) < 0.5*xw && fabs(p.y) < 0.5*yw) return flux/(xw*yw);
  else return 0.;  // do not use this function for fillXGrid()!
}

double 
SBBox::sinc(const double u) const {
  if (u<0.001 && u>-0.001)
    return 1.-u*u/6.;
  else
    return sin(u)/u;
}

DComplex
SBBox::kValue(Position<double> p) const{
  DComplex kval( sinc(0.5*p.x*xw)*sinc(0.5*p.y*yw), 0.);
  return kval*flux;
}

// Override fillXGrid so we can partially fill pixels at edge of box.
void 
SBBox::fillXGrid(XTable &xt) const {
  int N = xt.getN();
  double dx = xt.getDx(); // pixel grid size
  double norm = flux/xw/yw;
  
  // Pixel index where edge of box falls:
  int xedge = static_cast<int> ( ceil(xw / (2*dx) - 0.5) );
  int yedge = static_cast<int> ( ceil(yw / (2*dx) - 0.5) );
  // Fraction of edge pixel that is filled by box:
  double xfrac = xw / (2*dx) - xedge + 0.5;
  Assert(xfrac>0. && xfrac<=1.);
  double yfrac = yw / (2*dx) - yedge + 0.5;
  Assert(yfrac>0. && yfrac<=1.);
  if (xedge==0) xfrac = xw/dx;
  if (yedge==0) yfrac = yw/dx;

  double yfac;
  for (int iy = -N/2; iy < N/2; iy++) {
    if ( abs(iy) < yedge ) yfac = 0.;
    else if (abs(iy)==yedge) yfac = norm*yfrac;
    else yfac = norm;

    for (int ix = -N/2; ix < N/2; ix++) {
      if (yfac==0. || abs(ix)>xedge) xt.xSet(ix, iy ,0.);
      else if (abs(ix)==xedge) xt.xSet(ix, iy ,xfrac*yfac);
      else xt.xSet(ix,iy,yfac);
    }
  }
}

// Override x-domain writing so we can partially fill pixels at edge of box.
double
SBBox::fillXImage(Image<> I, double dx) const {
  double norm = flux/xw/yw;

  // Pixel index where edge of box falls:
  int xedge = static_cast<int> ( ceil(xw / (2*dx) - 0.5) );
  int yedge = static_cast<int> ( ceil(yw / (2*dx) - 0.5) );
  // Fraction of edge pixel that is filled by box:
  double xfrac = xw / (2*dx) - xedge + 0.5;
  Assert(xfrac>0. && xfrac<=1.);
  double yfrac = yw / (2*dx) - yedge + 0.5;
  Assert(yfrac>0. && yfrac<=1.);
  if (xedge==0) xfrac = xw/dx;
  if (yedge==0) yfrac = yw/dx;

  double totalflux = 0.;
  double xfac;
  for (int i = I.XMin(); i <= I.XMax(); i++) {
    if ( abs(i) > xedge ) xfac = 0.;
    else if (abs(i)==xedge) xfac = norm*xfrac;
    else xfac = norm;

    for (int j = I.YMin(); j <= I.YMax(); j++) {
      if (xfac==0. || abs(j)>yedge) I(i,j)=0.;
      else if (abs(j)==yedge) I(i,j)=xfac*yfrac;
      else I(i,j)=xfac;
      totalflux += I(i,j);
    }
  }
  // add image header to include "dx" info
  try {I.header()->erase("DK");} catch (ImageHeaderError &i) {}
  I.header()->replace("DX", dx);

  return totalflux * (dx*dx);
}



//
// SBLaguerre Class
//

// ??? Have not really investigated these:
double 
SBLaguerre::maxK() const {
  // Start with value for plain old Gaussian:
  double m=MAX(4., sqrt(-2.*log(ALIAS_THRESHOLD))) / sigma;
  // Grow as sqrt of order
  if (bvec.getOrder()>1) m *= sqrt(bvec.getOrder()/1.);
  return m;
}
double 
SBLaguerre::stepK() const {
  // Start with value for plain old Gaussian:
  double m= PI/MAX(4., sqrt(-2.*log(ALIAS_THRESHOLD))) / sigma;
  // Shrink as sqrt of order
  if (bvec.getOrder()>1) m /= sqrt(bvec.getOrder()/1.);
  return m;
 }

double
SBLaguerre::xValue(Position<double> p) const {
  LVector psi(bvec.getOrder());
  psi.fillBasis(p.x/sigma, p.y/sigma, sigma);
  double xval = bvec.dot(psi);
  return xval;
}

DComplex
SBLaguerre::kValue(Position<double> k) const {
  int N=bvec.getOrder();
  LVector psi(N);
  psi.fillBasis(k.x*sigma, k.y*sigma);  // Fourier[Psi_pq] is unitless
  // rotate kvalues of Psi with i^(p+q)
  // dotting b_pq with psi in k-space:
  double rr=0.;
  double ii=0.;
  {
    for (PQIndex pq(0,0); !pq.pastOrder(N); pq.nextDistinct()) {
      int j = pq.rIndex();
      double x = bvec[j]*psi[j] + (pq.isReal() ? 0 : bvec[j+1]*psi[j+1]);
      switch (pq.N() % 4) {
      case 0: 
	rr += x;
	break;
      case 1: 
	ii -= x;
	break;
      case 2: 
	rr -= x;
	break;
      case 3: 
	ii += x;
	break;
      }
    }  
  }
  // difference in Fourier convention with FFTW ???
  return DComplex(2*PI*rr, 2*PI*ii);
}

double
SBLaguerre::getFlux() const {
  double flux=0.;
  for (PQIndex pp(0,0); !pp.pastOrder(bvec.getOrder()); pp.incN())
    flux += bvec[pp].real();  // bvec[pp] is real, but need type conv.
  return flux;
}

void
SBLaguerre::setFlux(double flux_) {
  double newflux=flux_;
  if (getFlux()!=0) newflux /= getFlux();
  bvec.rVector() *= newflux;
  return;
}

// SBSersic Class 
// First need to define the static member that holds info on all the Sersic n's
SBSersic::InfoBarn SBSersic::nmap;

double
SBSersic::SersicInfo::kValue(double ksq) const {
  if (ksq<0.) 
    throw SBError("Negative k-squared passed to SersicInfo");
  if (ksq==0.) 
    return 1.;

  double lk=0.5*log(ksq);	// Lookup table is logarithmic

  if (lk<logkMin)
    return 1 + ksq*(kderiv2 + ksq*kderiv4);	// Use quartic approx at low k
  if (lk>=logkMax)
    return 0.;		// truncate the Fourier transform

  // simple linear interpolation to this value
  double fstep = (lk-logkMin)/logkStep;
  int index = static_cast<int> (floor(fstep));
  Assert(index < lookup.size()-1);
  fstep -= index;
  return lookup[index]*(1.-fstep) + fstep*lookup[index+1];
}

// Integrand class for the Hankel transform of Sersic
class SersicIntegrand {
public:
  SersicIntegrand(double n, double b_, double k_):
    invn(1./n), b(b_), k(k_) {}
  double operator()(double r) const {
    return r*exp(-b*pow(r, invn))*j0(k*r);
  }
private:
  double invn;
  double b;
  double k;
};


// Constructor to initialize Sersic constants and k lookup table
SBSersic::SersicInfo::SersicInfo(double n): inv2n(1./(2.*n)) {
  // Going to constraint range of allowed n to those I have looked at
  if (n<0.5 || n>4.2) throw SBError("Requested Sersic index out of range");

  // Formula for b from Ciotti & Bertin (1999)
  b = 2*n - (1./3.)
    + (4./405.)/n
    + (46./25515.)/(n*n)
    + (131./1148175.)/(n*n*n)
    - (2194697./30690717750.)/(n*n*n*n);

  double b2n = pow(b,2*n);  // used frequently here
  // The normalization factor to give unity flux integral:
  norm = b2n / (2*PI*n*tgamma(2.*n));

  // The quadratic term of small-k expansion:
  kderiv2 = -tgamma(4*n) / (4*b2n*tgamma(2*n)) ; 
  // And a quartic term:
  kderiv4 = tgamma(6*n) / (64*b2n*b2n*tgamma(2*n));

  /**cerr << "Building for n=" << n
	   << " b= " << b
	   << " norm= " << norm
	   << endl;
	   cerr << "Deriv terms: " << kderiv2 << " " << kderiv4 << endl; **/

  // When is it safe to use low-k approximation?  See when
  // quartic term is at threshold
  double lookupMin = 0.05;	// Default lower limit for lookup table
  const double kAccuracy=0.001;	// What errors in FT coefficients are acceptable?
  double smallK = pow(kAccuracy / kderiv4, 0.25);
  if (smallK < lookupMin) lookupMin = smallK;
  logkMin = log(lookupMin);

  // How far should nominal profile extend?
  // Estimate number of effective radii needed to enclose

  double xMax = 5.;	// Go to at least 5r_e
  {
    // Successive approximation method:
    double a=2*n;
    double z=a;
    double oldz=0.;
    int niter=0;
    const int MAXIT = 15;
    while ( abs(oldz-z)>0.01 && niter<MAXIT) {
      niter++;
      oldz = z;
      z = a - log(ALIAS_THRESHOLD*sqrt(2*PI*a)*(1+1./(12*a)+1./(288*a*a)))
	+(a-1)*log(z/a) + log(1 + (a-1)/z + (a-1)*(a-2)/(z*z));
    }
    double r=pow(z/b, n);
    if (r>xMax) xMax = r;
  }
  stepK = PI / xMax;

  // Going to calculate another outer radius for the integration of the 
  // Hankel transforms:
  double integrateMax=xMax;
  const double integrationLoss=0.001;
  {
    // Successive approximation method:
    double a=2*n;
    double z=a;
    double oldz=0.;
    int niter=0;
    const int MAXIT = 15;
    while ( abs(oldz-z)>0.01 && niter<MAXIT) {
      niter++;
      oldz = z;
      z = a - log(integrationLoss*sqrt(2*PI*a)*(1+1./(12*a)+1./(288*a*a)))
	+(a-1)*log(z/a) + log(1 + (a-1)/z + (a-1)*(a-2)/(z*z));
    }
    double r=pow(z/b, n);
    //**/cerr << "99.9% radius " << r <<endl;
    if (r>integrateMax) integrateMax = r;    
  }

  // Normalization for integral at k=0:
  double norm;
  const double INTEGRATION_RELTOL=0.0001;
  const double INTEGRATION_ABSTOL=1e-5;
  {
    SersicIntegrand I(n, b, 0.);
    // Integrate with at least 2^10 steps and up to 2^16:
    norm = Simp1d(I, 0., integrateMax, 
		  INTEGRATION_RELTOL, INTEGRATION_ABSTOL,
		  16, 10);
  }

  // Now start building the lookup table for FT of the profile.
  // Keep track of where the FT drops below ALIAS_THRESHOLD - this
  // will be our maxK.
  // Then extend the table another order of magnitude either in k
  //  or in FT, whichever comes first.
  logkStep = 0.05;
  // Here is preset range of acceptable maxK:
  const double MINMAXK = 10.;
  const double MAXMAXK = 50.; 
  maxK = MINMAXK;
  double lastVal=1.;
  double lk = logkMin;
  while (lk < log(maxK*10.) && lastVal>ALIAS_THRESHOLD/10.) {
    SersicIntegrand I(n, b, exp(lk));
    // Need to make sure we are resolving oscillations in the integral:
    double ncycles = integrateMax * exp(lk) / 2. / PI;
    int minSplit = static_cast<int> ( ceil(log(ncycles*4.) / log(2.) ) );
    minSplit = MAX(10, minSplit);
    //**/cerr << "minSplit " << minSplit << endl;
    double val = Simp1d(I, 0., integrateMax, 
			INTEGRATION_RELTOL, INTEGRATION_ABSTOL*norm,
			minSplit+6, minSplit);
    //**/cerr << "Integrate k " << exp(lk) << " result " << val/norm << endl;
    val /= norm;
    lookup.push_back(val);
    if (val >= ALIAS_THRESHOLD) maxK = MAX(maxK, exp(lk));
    logkMax = lk;
    lk += logkStep;
  }
  maxK = MIN(MAXMAXK, maxK);	// largest acceptable
}

// Integrand class for the flux integrals of Moffat
class MoffatFluxInt {
public:
  MoffatFluxInt(double beta_): beta(beta_) {}
  double operator()(double r) const {
    return r*pow(1.+r*r,-beta);
  }
private:
  double beta;
};

class MoffatFlux {
public:
  MoffatFlux(double beta): mfi(beta), target(0.) {}
  void setTarget(double target_) {target=target_;}
  double operator()(double r) const {
    return 2.*PI*Simp1d(mfi, 0., r) - target;
  }
private:
  MoffatFluxInt mfi;
  double target;
};


SBMoffat::SBMoffat(double beta_, double truncationFWHM,
		   double flux_, double re): beta(beta_),
					   flux(flux_),
					   ft(GTable<>::spline),
					   norm(1.),
					   rD(1.)
{
  //First, relation between FWHM and rD:
  FWHMrD = 2.* sqrt(pow(2., 1./beta)-1.);
  maxRrD = FWHMrD * truncationFWHM;

  {
    // Get flux and half-light radius
    MoffatFlux mf(beta);
    double fluxFactor = mf(maxRrD);
    solve::Solve<MoffatFlux> s(mf, 0.1, 2.);
    mf.setTarget(0.5*fluxFactor);
    double rerD = s.root();
    rD = re / rerD;
    norm = 1./fluxFactor;
  }

  // Make FFT's periodic at 4x truncation radius or 1.5x diam at ALIAS_THRESHOLD,
  // whichever is smaller
  stepKrD = 2*PI / MIN(4*maxRrD, 3.*sqrt(pow(ALIAS_THRESHOLD, -1./beta)-1));
  // And be sure to get at least 16 pts across FWHM when drawing:
  maxKrD = 16*PI / FWHMrD;

  // Get FFT by doing 2k transform over 2x the truncation radius
  // ??? need to do better here
  // ??? also install quadratic behavior near k=0?
  {
    const int N=2048;
    double dx = MAX(4*maxRrD, 64.) / N;
    XTable xt(N, dx);
    dx = xt.getDx();
    for (int iy=-N/2; iy<N/2; iy++)
      for (int ix=-N/2; ix<N/2; ix++) {
	double rsq = dx*dx*(ix*ix+iy*iy);
	xt.xSet(ix, iy, rsq<=maxRrD*maxRrD ? pow(1+rsq,-beta) : 0.);
      }
    KTable* kt = xt.transform();
    double dk = kt->getDk();
    double nn = kt->kval(0,0).real();
    for (int i=0; i<=N/2; i++) {
      ft.addEntry( i*dk, kt->kval(0,-i).real() / nn);
    }
    delete kt;
  }
}

DComplex
SBMoffat::kValue(Position<double> k) const {
  double kk = hypot(k.x, k.y)*rD;
  if (kk > ft.argMax()) return 0.;
  else return flux*ft(kk);
}
