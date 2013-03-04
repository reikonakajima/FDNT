// $Id: GLSimple2.cpp,v 1.9 2012/02/23 22:28:01 garyb Exp $
#include "GLSimple.h"
using namespace laguerre;

// Pixel management / masking routines.

// Class used to map (roughly) from world coordinates to pixels
void 
CrudeMap::pixToWorld(double px, double py, double& wx, double& wy) {
    px -= p0[0]; py -= p0[1];
    wx = dWdP(0,0)*px + dWdP(0,1)*py + w0[0];
    wy = dWdP(1,0)*px + dWdP(1,1)*py + w0[1];
  }
void 
CrudeMap::worldToPix(double wx, double wy, double& px, double& py) {
    wx -= w0[0]; wy -= w0[1];
    px = dPdW(0,0)*wx + dPdW(0,1)*wy + p0[0];
    py = dPdW(1,0)*wx + dPdW(1,1)*wy + p0[1];
  }
Ellipse
CrudeMap::worldToPix(const Ellipse& ew) {
  Ellipse ep = Ellipse::fromMatrix(dPdW * ew.getMatrix());
  Position<double> wx0 = ew.getX0();
  double px0, py0;
  worldToPix(wx0.x, wx0.y, px0, py0);
  ep.setX0( Position<double>(px0,py0));
  return ep;
}
// Initialize by getting local map from point in an image close to
// the specified world coord system.  Use near edge if off image.
CrudeMap::CrudeMap(const Image<double> xw, const Image<double> yw,
				double xw0, double yw0): 
  dWdP(2,2), dPdW(2,2), w0(2), p0(2) {
  Assert (xw.getBounds()==yw.getBounds());
  // Start with pixel at center of image
  int ipx, ipy;
  {
    Position<int> ip = xw.getBounds().center();
    ipx = ip.x; ipy = ip.y;
  }
  int iter=0;
  const int MAX_ITER=6;	// allowed iterations to localize to 1 pixel	
  do {
    iter++;
    if (iter > MAX_ITER)
      throw GLError("Too many iterations in CrudeMap");
    p0[0] = ipx; p0[1] = ipy;
    w0[0] = xw(ipx,ipy);
    w0[1] = yw(ipx,ipy);
    // x derivs
    if (ipx==xw.XMin()) {
      // Right-side derivs:
      dWdP(0,0) = xw(ipx+1,ipy)-xw(ipx,ipy);
      dWdP(1,0) = yw(ipx+1,ipy)-yw(ipx,ipy);
    } else if (ipx==xw.XMax()) {
      // Left-side derivs:
      dWdP(0,0) = xw(ipx,ipy)-xw(ipx-1,ipy);
      dWdP(1,0) = yw(ipx,ipy)-yw(ipx-1,ipy);
    } else {
      // symmetric derivs:
      dWdP(0,0) = 0.5*(xw(ipx+1,ipy)-xw(ipx-1,ipy));
      dWdP(1,0) = 0.5*(yw(ipx+1,ipy)-yw(ipx-1,ipy));
    }

    // y derivs
    if (ipy==xw.YMin()) {
      // Right-side derivs:
      dWdP(0,1) = xw(ipx,ipy+1)-xw(ipx,ipy);
      dWdP(1,1) = yw(ipx,ipy+1)-yw(ipx,ipy);
    } else if (ipy==xw.YMax()) {
      // Left-side derivs:
      dWdP(0,1) = xw(ipx,ipy)-xw(ipx,ipy-1);
      dWdP(1,1) = yw(ipx,ipy)-yw(ipx,ipy-1);
    } else {
      // symmetric derivs:
      dWdP(0,1) = 0.5*(xw(ipx,ipy+1)-xw(ipx,ipy-1));
      dWdP(1,1) = 0.5*(yw(ipx,ipy+1)-yw(ipx,ipy-1));
    }

    // Save inverse
    dPdW = dWdP.inverse();

    // Map desired center to a pixel
    double px, py;
    worldToPix(xw0, yw0, px, py);
      
    int ipxnew = static_cast<int> (floor(px+0.5));
    int ipynew = static_cast<int> (floor(py+0.5));
    ipxnew = MIN(ipxnew, xw.XMax());
    ipxnew = MAX(ipxnew, xw.XMin());
    ipynew = MIN(ipynew, xw.YMax());
    ipynew = MAX(ipynew, xw.YMin());
    if (abs(ipxnew-ipx)<=1 && abs(ipynew-ipy)<=1 ) return;
    ipx = ipxnew;
    ipy = ipynew;
  } while (true);
}

// First routine constructs mask and obtains all data within mask.
// 
template <class T>
bool
GLSimple<T>::acquireData() {
  // Skip if we're already set:
  if (dataAreValid) 
    return (!(flags & (OutOfBounds + TooLarge)));

  // else clear flags and roll:
  flags &= ~(Edge + OutOfBounds + TooLarge);

  // Select ellipse for mask:
  maskBasis = basis;
  double useSigma = maskSigma;
  if (useSigma <= 0.) {
        useSigma = MAX(minimumMaskSigma, 
		       sqrtMaskFactor*sqrt(static_cast<double> (order)));
  }
  maskBasis.setMu( maskBasis.getMu() + log(useSigma) );
  if (maskBasis.getMu() > maxMaskMu) {
    flags |= TooLarge;
    dataAreValid = false;
    return false;
  }

  // Initialize vectors for data points
  vector<double> xx;
  vector<double> yy;
  vector<double> dd;
  vector<double> ss;

  // OutOfBounds set until we find 1 image in bounds
  flags |= OutOfBounds;

  // Loop through images:

  bool buildMaps = maps.empty();

  int expno=0;
  double xw0 = basis.getX0().x;
  double yw0 = basis.getX0().y;
  for ( typename list<FitExposure<T> >::const_iterator iexp = exposures.begin();
	iexp != exposures.end();
	++iexp, ++expno) {
    if (buildMaps) maps.push_back(CrudeMap( iexp->xworld, iexp->yworld, xw0, yw0));

    Bounds<int> coverage = iexp->sci.getBounds() & iexp->wt.getBounds();

    dbg << "coverage: " << coverage << endl;

    // If center is on the image && valid data, clear OutOfBounds
    double xp0, yp0;
    maps[expno].worldToPix(xw0, yw0, xp0, yp0);
    int ipx0 = static_cast<int> (floor(xp0+0.5));
    int ipy0 = static_cast<int> (floor(yp0+0.5));
    if ( coverage.includes(ipx0, ipy0) && iexp->wt(ipx0,ipy0)>0.)
      flags &= ~OutOfBounds;

    // Map the ellipse into pixel coordinates, roughly
    Ellipse pixBasis = maps[expno].worldToPix(maskBasis);
    // Find pixel-coord Bounds that will enclose the maskBounds
    Bounds<double> db = pixBasis.range();
    Bounds<int> pixBounds( static_cast<int> (floor(db.getXMin())),
			   static_cast<int> (ceil(db.getXMax())),
			   static_cast<int> (floor(db.getYMin())),
			   static_cast<int> (ceil(db.getYMax())));
    dbg  << expno << " pixBasis: " << pixBasis
	 << " range " << pixBounds
	 << " img bounds " << coverage
	 << endl;

    // Set Edge flag if pixels outside image would contribute
    if ( !coverage.includes(pixBounds)) 
      flags |= Edge;

    pixBounds = pixBounds & coverage;

    // Search through pixelBounds for valid data
    if (pixBounds) 
      for (int iy = pixBounds.getYMin(); iy<=pixBounds.getYMax(); iy++)
	for (int ix = pixBounds.getXMin(); ix<=pixBounds.getXMax(); ix++)
	  if (iexp->wt(ix,iy) > 0.) {
	    double invs=sqrt(iexp->wt(ix,iy));
	    xx.push_back(iexp->xworld(ix,iy));
	    yy.push_back(iexp->yworld(ix,iy));
	    dd.push_back((iexp->sci(ix,iy)-iexp->sky) * invs);
	    ss.push_back( invs / iexp->sbScaleFactor);
	  }
  }
  // Copy data into DVectors
  const int npts = xx.size();
  if (x) delete x;
  x = new DVector(npts);
  std::copy(xx.begin(), xx.end(), x->begin());
  if (y) delete y;
  y = new DVector(npts);
  std::copy(yy.begin(), yy.end(), y->begin());
  if (dsig) delete dsig;
  dsig = new DVector(npts);
  std::copy(dd.begin(), dd.end(), dsig->begin());
  if (invsig) delete invsig;
  invsig = new DVector(npts);
  std::copy(ss.begin(), ss.end(), invsig->begin());

  dataAreValid = true;
  return (!(flags & (OutOfBounds + TooLarge)));
}

template <class T>
void
GLSimple<T>::maskRatios(double& smallest, double& largest) {
  // Select ellipse we'd want now:
  Ellipse tryBasis = basis;
  double useSigma = maskSigma;
  if (useSigma <= 0.) {
        useSigma = MAX(minimumMaskSigma, 
		       sqrtMaskFactor*sqrt(static_cast<double> (order)));
  }
  tryBasis.setMu( tryBasis.getMu() + log(useSigma) );

  // Tramp around edge of new mask, ask what it's "radius" is
  // in the old mask:

  const int nSteps=16;
  smallest = 100.;
  largest = 0.;
  for (int istep=0; istep<nSteps; istep++) {
    Position<double> xin( cos(istep*2*PI/nSteps),
			  sin(istep*2*PI/nSteps));
    Position<double> xout = maskBasis.inv(tryBasis.fwd(xin));
    double r = hypot(xout.x, xout.y);
    smallest = MIN(smallest, r);
    largest = MAX(largest, r);
  }
}

