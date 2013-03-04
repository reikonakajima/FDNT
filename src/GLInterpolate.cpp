// $Id$
// Test code for interpolating locations of missing data using GL fit to fill gaps.
#include "SBParse.h"
#include "StringStuff.h"
#include "Std.h"
#include "FITSImage.h"
#include "Random.h"
#include "Laguerre.h"
#include "GLSimple.h"
#include <list>

using namespace sbp;
using namespace img;
using namespace laguerre;

string usage = "GLInterpolate: demonstrate interpolation over bad pixels using GL fitting\n"
               "usage:  GLInterpolate <S/N>  <order> <radius>\n"
               " <S/N>  is signal-to-noise ratio to draw image (with matched filter)\n"
               "<order> is order of Gauss-Laguerre fit\n"
               "<radius> is initial guess for size of GL basis.  Initial guess is centered\n"
               "        at origin with zero ellipticity.\n"
               " stdin has first line describing the shape to draw using SBParse syntax\n"
               "       and succeeding lines are pairs of integers giving coords of bad pixels\n"
               " FITS images stamp.fits and wt.fits will be produced holding interpolated data\n"
               " and weight map, respectively.";
int
main(int argc, char *argv[])
{
  if (argc!=4) {
    cerr << usage << endl;
    exit(1);
  }
  double snr = atof(argv[1]);
  int order = atof(argv[2]);
  double startRadius = atof(argv[3]);
  try {
    // Create image using SBParse
    string buffer;
    cerr << "Enter SB specification for parser: " << endl;
    stringstuff::getlineNoComment(cin, buffer);
    SBProfile* sbp = SBParse(buffer);
    //**/cerr << "Back from parse " << endl;
    Bounds<int> b(-32,31,-32,31);
    Image<> stamp(b);
    sbp->draw(stamp,1.);
    //**/cerr << "Back from draw " << endl;
    b = stamp.getBounds();
    // The stamp image now contains a noiseless version of the galaxy
    double fluxPerfect = sbp->getFlux();
    delete sbp;

    /**/cerr << "stamp bounds " << stamp.getBounds() << endl;

    // Choose a noise level to give desired matched-filter S/N
    double squaresum = 0;
    for (int iy = b.getYMin(); iy<=b.getYMax(); iy++)
      for (int ix = b.getXMin(); ix<=b.getXMax(); ix++)
	squaresum+=stamp(ix,iy)*stamp(ix,iy);

    double noise = sqrt(squaresum)/snr;
    // This value of Gaussian noise will satisify (sum (signal/noise)^2) = snr^2,
    // where the sum is over all pixels.  Now add a realization of Gaussian noise to
    // the image and make a (uniform) weight map from inverse variance.
    {
      ran::GaussianDeviate g;
      for (int iy = b.getYMin(); iy<=b.getYMax(); iy++)
	for (int ix = b.getXMin(); ix<=b.getXMax(); ix++)
	  stamp(ix,iy) = stamp(ix,iy) + g*noise;
    }
    Image<> wt(b, pow(noise, -2.));

    // punch holes in weight map
    vector<int> xBad;
    vector<int> yBad;
    cerr << "Enter bad pixel coords, one pair per line, end with EOF:" << endl;
    while (getlineNoComment(cin, buffer)) {
      int ix, iy;
      istringstream iss(buffer);
      if (!(iss >> ix >> iy)) {
	cerr << "Bad specification for coordinates of defect: "
	     << buffer
	     << endl;
	exit(1);
      }
      if (b.includes(ix,iy)) {
	wt(ix, iy) = 0.;
	xBad.push_back(ix);
	yBad.push_back(iy);
      }
    }

    // Produce (trivial) coordinate map for the postage stamp
    Image<double> xw(b);
    Image<double> yw(b);
    for (int iy = b.getYMin(); iy<=b.getYMax(); iy++)
      for (int ix = b.getXMin(); ix<=b.getXMax(); ix++) {
	xw(ix,iy) = ix;
	yw(ix,iy) = iy;
    }

    // Execute the GL fitting:
    FitExposure<> fe(stamp, wt, xw, yw);
    Ellipse basis(0., 0., log(startRadius)); // starting basis ellipse for GL fitting
    GLSimple<> gl(fe, basis, order);
    gl.solve();
    cout << "Flags: " << gl.getFlags() << endl;
    basis = gl.getBasis();
    cout << "Basis: " << basis << endl;
    LVector bvec = gl.getB();
    double fluxModel = bvec.flux();

    // Calculate missing-pixel values
    double missingFlux=0.;
    double scaleFactor = exp(basis.getMu());

    for (int ipix = 0; ipix<xBad.size(); ipix++) {
      int ix = xBad[ipix];
      int iy = yBad[ipix];
      // First map the pixel coordinates so that basis ellipse becomes unit circle:
      Position<double> xunit = basis.inv(Position<double>(xw(ix,iy),yw(ix,iy)));
      LVector psi(bvec.getOrder());
      // The fillBasis method sets the LVector elements to the size of their corresponding
      // Gauss-Laguerre basis functions at the position requested.  The scaleFactor, which is
      // the determinant of the Ellipse transformation, is needed here for proper normalization
      // of the basis functions.
      psi.fillBasis(xunit.x, xunit.y, scaleFactor);
      /**/cerr << "At (" << ix << "," << iy << ") Old " << stamp(ix,iy);
      stamp(ix,iy) = bvec.dot(psi);
      /**/cerr << " interpolated " << stamp(ix,iy) << endl;
      missingFlux += stamp(ix,iy);
    }
    // Report fraction of flux interpolated:
    cout << "True flux " << fluxPerfect
	 << " model flux " << fluxModel
	 << " interpolated fraction " << missingFlux / fluxModel
	 << endl;

    // Draw resultant images to disk.
    stamp.shift(1,1);
    FITSImage<>::writeToFITS("stamp.fits", stamp);
    wt.shift(1,1);
    FITSImage<>::writeToFITS("wt.fits", wt);

  } catch (std::runtime_error& m) {
    quit(m,1);
  }
  exit(0);
}

