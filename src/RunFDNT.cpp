// RunFDNT.cpp
// Run FDNT on images with given galaxy and PSF images

#include "RunFDNT.h"
#include "SBPixel.h"

using namespace laguerre;
using namespace sbp;
using namespace astrometry;


namespace fdnt {

char const* greet()
{
  return "hello, world";
}

template <typename T>
FDNTShapeData RunFDNT(const Image<T>& gal_image, const Image<T>& psf_image,
		      const Image<T>& weight_image,
		      double x_pix, double y_pix,
		      double a_wc, double b_wc, double pa_wc,
		      double r_pix,    // FLUX_RADIUS in pixels, not wcs
		      double ee50psf,  // PSF half-light radius
		      double bg, int order, double sky) {

  // necessary input:
  // sky: sky_level, order: FDNT GL-fit order (0 for FFT)
  double maskSigma = -1.0;  // GL and FDNT mask sigma (-1 auto)
  string weightScaleKey = "WTSCALE"; // Scaling factor for weight map keyword in WCS header
  // This is the scale that transforms weight proportional to 1/sig^2 to weights EQUAL to
  // 1/sig^2 of the respective single frame.
  // Note that single frame rescaling is happening on top of this on both the science and
  // weight frames
  string fluxScaleKey = "FLXSCALE";
  int psfOrder;  // maximum order of polynomial psf variation used; -1 for order of PSFEx model
  int stampSize = 64;  // Pixel size of img postage stamps

  // GL interpolation of missing values
  double maxBadPixels = 0.1;  // Maximum fraction of bad pixels in postage stamp
  int interpolationOrder = 4; // GL order for bad pixel interpolation
  double maxBadFlux = 0.05;   // Maximum fraction of model flux in bad pixels

  // g1Col: g1 shape estimate (use native shape if 0)
  // g2Col: g2 shape estimate (use native shape if 0)

  // output:
  FDNTShapeData results;

  try {

    //PSFExModel *model;
    Image<T> sci = gal_image.duplicate();
    /*  /// something to include in the future
	img::ImageHeader h;
	ifstream mapfs(wcsName.c_str());
	if (!mapfs) {
	    cerr << "Could not open astrometry file " << wcsName
		 << "; trying FITS header in science frame" << endl;
	    h = *(sci.header());
	}
	else {
	    h = img::HeaderFromStream(mapfs);
	}
	SCAMPMap *fullmap = new SCAMPMap(h,0);
	double xw, yw;
	fullmap->toWorld(0,0,xw,yw);
        cerr << "# WCS " << xw << " " << yw << " is the (0,0) pixel WC w.r.t. the TP" << endl;
    */

    Image<T> wt = weight_image.duplicate();
    double weightScale = 1.0;  // SCAMP weight scale
    double fluxScale = 1.0;

    Bounds<int> chip_bounds = sci.getBounds();
    Bounds<int> safe_chip_bounds = chip_bounds;
    safe_chip_bounds.addBorder(-2);

    // the inverse variance scales with flux in the following way:
    weightScale = weightScale / (fluxScale*fluxScale);
    for (int iy=chip_bounds.getYMin(); iy<=chip_bounds.getYMax(); iy++)
      for (int ix=chip_bounds.getXMin(); ix<=chip_bounds.getXMax(); ix++) {
	// weight map just rescaled sky inverse-variance:
	wt(ix,iy)  = wt(ix,iy) * weightScale;
	sci(ix,iy) = sci(ix,iy) * fluxScale;
      }

    UnweightedShearEstimator se;
    //tmv::SymMatrix<double> covE(2);

    // DEBUG: REMOVE
    cerr << "# x_pix y_pix eta1 eta2 sig1 sig2 cov12 mu egFix fdFlags considered_success"
	 << endl;

    double x_wc = x_pix, y_wc = y_pix;

    /*  /// future project
	double x_wc, y_wc; // tangent plane coordinates
	// deproject spherical coordinates to xi (x_wc) and eta (y_wc)
	SphericalICRS point(ra0*DEGREE, dec0*DEGREE);
	TangentPlane tp(point, fullmap->projection());
	tp.getLonLat(x_wc, y_wc);
	x_wc /= DEGREE;
	y_wc /= DEGREE;
	try {
	  //cerr << "# object at WCS " << x_wc << " " << y_wc << endl; // DEBUG
	  fullmap->toPix(x_wc, y_wc, x_pix, y_pix);
	}
	catch  (AstrometryError& e) {
	  cerr << e.what() << endl;
	  cerr << "(b) processing object at (RA,dec)=(" << ra0 << "," << dec0 << ")" << endl;
	  cerr << "toPix failure" << endl;
	  exit(0);
	}
    */
    if (!safe_chip_bounds.includes(int(x_pix), int(y_pix)))
      throw MyException("x_pix, y_pix not within image bounds");

    /*  /// future project
	Matrix22 J = fullmap->dWorlddPix(x_pix,y_pix);
	// distorted approximate pixel scale: how much the scales are enlarged by WC transform
	double rescale = sqrt(fabs(J(0,0)*J(1,1)-J(0,1)*J(1,0)));
    */
    double rescale = 1.0;
    Ellipse psfBasis(0., 0., log(ee50psf/1.17));

    const Quintic quintic(1e-4);
    InterpolantXY quintic_2d(quintic);

    SBPixel psfWCS(psf_image, quintic_2d);
    /* // generate PSF image in WC  /// future project
        SBDistort psfWCS(*(model->sb()),J(0,0),J(0,1),J(1,0),J(1,1));
	cerr << "dWorld/dPix " << J;  // J comes with its own endl
	cerr << ee50psf << endl;
    */
    // sets flux of basis correctly to have flux normalization after distortion
    psfWCS.setFlux(1.);
    double dx = ee50psf / 2.35;  // for proper sampling
    Image<T> ipsf = psfWCS.draw(dx);  // future project: generate from model

    // Measure PSF GL size & significance
    Image<T> psfwt(ipsf.getBounds());
    psfwt = pow(1e-5/dx, -2.);
    psfBasis.setMu( psfBasis.getMu() - log(dx));

    GLSimple<> gl(ipsf, psfwt, psfBasis, 4);
    bool psf_success = gl.solve();
    if (!psf_success) {
      std::ostringstream oss;
      oss << "# Failed measuring psf, flags " << gl.getFlags()
	  << "; ignoring object" << endl;
      throw MyException(oss.str().c_str());
    }
    psfBasis = gl.getBasis();
    psfBasis.setMu( psfBasis.getMu() + log(dx));

    // save measurement results on observed galaxy properties
    results.psf_flags = gl.getFlags();
    if (psf_success) {
      results.psf_e1 = psfBasis.getS().getE1();
      results.psf_e2 = psfBasis.getS().getE2();
      results.psf_sigma = exp(psfBasis.getMu());
      gl.b00(results.psf_b00, results.psf_b00_var);
      results.psf_order = gl.getOrder();
      if (results.psf_order >= 4) {
	results.psf_b22 = gl.getB()(2,2).real();
      }
      results.psf_chisq = gl.getChisq();
      results.psf_DOF = gl.getDOF();
    }

    // setup PSF information for deconvolution
    PSFInformation psfinfo(psfWCS, psfBasis);

    double e1start = (a_wc*a_wc-b_wc*b_wc)/(a_wc*a_wc+b_wc*b_wc);
    double e2start = e1start * sin(2*pa_wc*PI/180.);
    e1start *= cos(2*pa_wc*PI/180.);
    // all in wcs pixel coordinates
    Ellipse sexE(e1start, e2start, log(r_pix*rescale), x_wc, y_wc);

    Shear initShear;
    Ellipse initE;
    double g1_wc;
    double g2_wc;

    /*  /// future project: get external initial shape estimate
	if (g1Col > 0 && g2Col > 0) {
	  // read in g1, g2 (if available)
	  g1_wc = atof(readvals[g1Col-1].c_str());
	  g2_wc = atof(readvals[g2Col-1].c_str());
	  double gabs = sqrt(g1_wc*g1_wc+g2_wc*g2_wc);
	  if (gabs > 0.95) {  // sanitize shear
	    g1_wc=g1_wc*0.95/gabs;
	    g2_wc=g2_wc*0.95/gabs;
	  }
	  initShear.setG1G2(g1_wc, g2_wc);
	  initShear.getE1E2(e1start, e2start);
	  initE = Ellipse(e1start, e2start, log(r_pix*rescale), x_wc, y_wc); // all in wcs
    */

    // approximate g1, g2 with native shape (if external g1, g2 not available)
    initE = sexE;

    FitExposure<T>* fep;

    // set-up bounds for postagestamp
    int xi0 = x_pix-stampSize/2;
    if (xi0 < chip_bounds.getXMin())
      xi0 = chip_bounds.getXMin();
    if ((xi0 + stampSize+1) > chip_bounds.getXMax())
      xi0 = chip_bounds.getXMax() - stampSize - 1;
    int yi0 = y_pix-stampSize / 2;
    if (yi0 < chip_bounds.getYMin())
      yi0 = chip_bounds.getYMin();
    if ((yi0 + stampSize+1) > chip_bounds.getYMax())
      yi0 = chip_bounds.getYMax() - stampSize - 1;

    Bounds<int> stamp(xi0,xi0+stampSize-1, yi0, yi0+stampSize-1);
    results.image_bounds = stamp;
    Image<float> scistamp = sci.subimage(stamp);
    Image<float> wtstamp  = wt.subimage(stamp);
    Image<double> xwstamp(stamp);
    Image<double> ywstamp(stamp);

    try {
      for (int iy=stamp.getYMin(); iy<=stamp.getYMax(); iy++)
	{
	  for (int ix=stamp.getXMin(); ix<=stamp.getXMax(); ix++) {
	    double dxw = ix, dyw = iy;  // pixel coordinate is world coordinate
	    //fullmap->toWorld(ix,iy,dxw,dyw);  /// future project
	    xwstamp(ix,iy) = dxw;
	    ywstamp(ix,iy) = dyw;
	    /*  /// deweight any irrelevant pixel (future project?)
		if (segstamp(ix,iy) && segstamp(ix,iy)!=segId) { //it's another object!
		  for (int iiy=max(stamp.getYMin(),iy-segmentationPadding);
		       iiy<=min(stamp.getYMax(),iy+segmentationPadding); iiy++)
		    for (int iix=max(stamp.getXMin(),ix-segmentationPadding);
			 iix<=min(stamp.getXMax(),ix+segmentationPadding); iix++)
		      wtstamp(iix,iiy)=0.;
	    */
	  }
	}
      fep = new FitExposure<>(scistamp,
			      wtstamp,
			      xwstamp,
			      ywstamp, 0, sky+bg);
      // stack background with sky-subtracted single frames ==
      // photometric local background in stack flux scale
    } catch (...) {
      throw MyException("# fail: could not get postage stamp");
    }

    double meanweight=0.;
    vector< Position<int> > bp = fep->badPixels(meanweight);
    if (bp.size()>maxBadPixels*stampSize*stampSize)
      throw MyException ("too many bad pixels");

    // has bad pixels, but not too many to begin with: do GL interpolation
    if (bp.size() > 0) {
      GLSimple<> gal(*fep, sexE, interpolationOrder);
      if (!gal.solve())
	throw MyException("GL fit failed");

      LVector bvec = gal.getB();
      double fluxModel = bvec.flux();
      Ellipse basis = gal.getBasis();
      double missingFlux = 0.;
      double scaleFactor = exp(basis.getMu());

      for (vector< Position<int> >::iterator it=bp.begin(); it<bp.end(); it++) {
	LVector psi(bvec.getOrder());
	Position<double> xunit =
	  basis.inv(Position<double>(xwstamp((*it).x,(*it).y),
				     ywstamp((*it).x,(*it).y)));
	psi.fillBasis(xunit.x, xunit.y, scaleFactor);
	scistamp((*it).x,(*it).y)=bvec.dot(psi);
	wtstamp((*it).x,(*it).y)=meanweight;
	// The interpolated pixel is assumed to have the mean weight of the good pixels;
	// this makes sense because in the Fourier code homogeneous uncertainties are
	// assumed
	missingFlux += scistamp((*it).x,(*it).y);
	scistamp((*it).x,(*it).y)+=sky+bg;
      }
      missingFlux *= rescale*rescale; // scale flux with coordinate system

      if (missingFlux/fluxModel > maxBadFlux)
	throw MyException("bad pixels have too high flux fraction");
    }

    FDNT<> fd(*fep, psfinfo, initE, order);
    fd.setMaskSigma(maskSigma);
    fd.GLAll();
    bool success = fd.prepare();

    Shear targetS = fd.getBasis().getS();
    double prob;
    tmv::SymMatrix<double> covE(2,2);
    if (success) {
      try {
	targetS = fd.shape2(prob, covE);  // THIS IS THE ACTUAL MEASUREMENT!!
      }
      catch (...) {
	  ;  // do nothing
      }
      // ??? Make flag mask an input parameter ???
      success = !(fd.getFlags() & (DidNotConverge + Singularity + OutOfBounds + TooLarge
				   + UnderSampled + TooElliptical + GLFailure));
    }

    double eta1, eta2;
    double g1, g2;
    targetS.getEta1Eta2(eta1, eta2);
    targetS.getG1G2(g1, g2);
    double egFix = 0.;
    double sig1 = 0., sig2 = 0.;
    double cov12 = 0.;

    if (success) {
      try {
	egFix = fd.shrinkResponse(targetS);
	sig1 = sqrt(covE(0,0));
	sig2 = sqrt(covE(1,1));
	cov12 = covE(0,1);
	se.add(targetS, egFix, sig1, sig2);
      }
      catch (...)
      {
	throw MyException("an error occured on previously successful measurement");
      }
    } else {
      cerr << "# fail: some obscure thing in FDNT::prepare() happens: "<< fd.getFlags() << endl;
    }

    double mu = fd.getBasis().getMu();

    // save measurement results on intrinsic galaxy properties
    results.intrinsic_flags = fd.getFlags();
    if (success) {
      results.intrinsic_e1 = fd.getBasis().getS().getE1();
      results.intrinsic_e2 = fd.getBasis().getS().getE2();
      results.intrinsic_e1_var = covE(0,0);
      results.intrinsic_e2_var = covE(1,1);
      results.intrinsic_e1e2_covar = covE(0,1);
      results.intrinsic_sigma = exp(mu);
      results.shrink_response = egFix;
      results.evaluation_count = fd.getEvaluationCount();
      results.e_trial_count = fd.getETrialCount();
    }

    // DEBUG: REMOVE
    cerr << x_pix
	 << " " << y_pix
	 << " " << eta1
	 << " " << eta2
	 << " " << sig1
	 << " " << sig2
	 << " " << cov12
	 << " " << mu
	 << " " << egFix
	 << " " << fd.getFlags()
	 << " " << success
	 << endl;

    delete fep;

    //delete model;   /// future project, to include psfmodel and astrometry
    //delete fullmap;

    // Print out mean shear estimate
    Shear S(se);
    S.getG1G2(g1,g2);
    se.sigmaE(sig1,sig2, false);

    // DEBUG: REMOVE
    cerr << fixed << setprecision(6);
    // Approximate the reduced-shear error as 1/2 of the distortion error
    cerr << "# Means: " << g1 << " +- " << sig1/2
	 << " " << g2 << " +- " << sig2/2
	 << endl;
  } catch (std::runtime_error &m) {
    cerr << m.what() << endl;
    quit(m,1);
  }
  return results;
}

// instantiate for expected types

template FDNTShapeData RunFDNT<float>(const Image<float>&, const Image<float>&, const Image<float>&,
				      double, double, double, double, double, double, double,
				      double, int, double);
/*
template FDNTShapeData RunFDNT<double>(const Image<double>&, const Image<double>&,
				       const Image<double>&,
				       double, double, double, double, double, double, double,
				       double, int, double);
*/

} // namespace fdnt
