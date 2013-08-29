// PSFEx_KiDS.cpp
// Generate PSFEx psf models and measure shape; generate catalog

#include "FDNT.h"
#include "StringStuff.h"
#include "Pset.h"
#include "FITSImage.h"
#include <fstream>
#include "ShearEstimator.h"
#include "EnclosedFluxRadius.h"
#include "PSFEx.h"
#include "Astrometry.h"
#include "Image.h"
#include "HeaderFromStream.h"
#include "SCAMPMap.h"

//#define DEBUGFDNTPSFEX
//#define CHECKPLOTS

using namespace laguerre;
using namespace sbp;
using namespace astrometry;

int
main(int argc, 
     char *argv[]) {
    
    // Parameters to set:
    // columns for necessary input:
    int idCol;
    int segIdCol;
    int raCol;
    int decCol;
    int magCol;
    int bgCol;
    int rCol;              // half-light radius (pix)
    int aCol, bCol, paCol; // WC observed ellipse
    int g1Col;
    int g2Col;
    int fwdColStart;
    int fwdColEnd;
    
    // The fit:
    int order;
    double maskSigma;
    
    // Input files
    string fitsName;
    string catName;
    string psfName;
    int    psfOrder;
    string weightName;
    string wcsName;
    string segmentationName;
    int    segmentationPadding;
    
    string weightScaleKey;
    string fluxScaleKey;
    
    // Data properties
    double sky;
    //double rn;
    //double gain;
    int stampSize;
    
    // GL interpolation of missing values
    double maxBadPixels;
    int interpolationOrder;
    double maxBadFlux;
    
    Pset parameters;
    {
	const int def=PsetMember::hasDefault;
	const int low=PsetMember::hasLowerBound;
	const int up=PsetMember::hasUpperBound;
	const int lowopen = low | PsetMember::openLowerBound;
	const int upopen = up | PsetMember::openUpperBound;
	
	parameters.addMemberNoValue("FILES:",0,
				    "Input data");
	parameters.addMember("fitsName", &fitsName, def,
			     "Image FITS file", "set0001.fit");
	parameters.addMember("catName", &catName, def,
			     "Input SExtractor catalog", "set0001.scat");
	parameters.addMember("psfName", &psfName, def,
			     "Input PSFEx file", "model.txt");
	parameters.addMember("psfOrder", &psfOrder, def | low,
			     (std::string("Maximum order of polynomial psf variation used; ")+
			      std::string("-1 for order of PSFEx model")).c_str(), -1, -1);
	parameters.addMember("weightName", &weightName, def, 
			     "Input weight map proportional to inverse variance", "weight.fit");
	parameters.addMember("wcsName", &wcsName, def,
			     "World coordinate header file", "wcs.txt");
	parameters.addMember("weightScaleKey", &weightScaleKey, def,
			     "Scaling factor for weight map (e.g. from Swarp) keyword in WCS header",
			     "WTSCALE");
	// This is the scale that transforms weight proportional to 1/sig^2 to weights EQUAL to 
	// 1/sig^2 of the respective single frame.
	// Note that single frame rescaling is happening on top of this on both the science and 
	// weight frames
	parameters.addMember("fluxScaleKey", &fluxScaleKey, def,
			     "Scaling factor for flux (e.g. from WCS header) keyword in WCS header",
			     "FLXSCALE");
	parameters.addMember("segmentationName", &segmentationName, def,
			     "segmentation map", "segment.fits");
	parameters.addMember("segmentationPadding", &segmentationPadding, def | low,
			     "padding around segmentation maps [pix]", 0, 0);
	parameters.addMember("stampSize", &stampSize, def | low,
			     "Pixel size of img postage stamps", 64, 3);
	
	parameters.addMemberNoValue("DATA PROPERTIES:",0,
				    "Input data");
	parameters.addMember("sky", &sky, def,
			     "sky level", 0.);
	//parameters.addMember("rn", &rn, def | low,
	//		   "Read noise, i.e. RMS at 0 ADU", 1000., 0.);
	//parameters.addMember("gain", &gain, def | low,
	//		   "Gain, e/ADU for Poissson, 0 for no Poisson", 1., 0.);
	
	parameters.addMemberNoValue("FDNT:",0,
				    "Fitting parameters:");
	parameters.addMember("order", &order, def,
			     "FDNT GL-fit order (0 for FFT)", 0);
	parameters.addMember("maskSigma", &maskSigma, def,
			     "GL and FDNT mask sigma (-1 auto)", -1. /*4.*/); // changed 21.11.12
	
	parameters.addMemberNoValue("CATALOG TRANSLATION:",0,
				    "Columns holding input data");
	parameters.addMember("idCol", &idCol, def | low,
			     "Object ID", 1, 1);
	parameters.addMember("segIdCol", &segIdCol, def | low,
			     "Object ID in segmentation map; determined automatically if 0", 1, 0);
	parameters.addMember("raCol", &raCol, def | low,
			     "RA centroid", 2, 1);
	parameters.addMember("decCol", &decCol, def | low,
			     "DEC centroid", 3, 1);
	parameters.addMember("magCol", &magCol, def | low,
			     "Magnitude", 4, 1);
	parameters.addMember("bgCol", &bgCol, def | low,
			     "photometric background flux (zero if bgCol==0)", 5, 0);
	parameters.addMember("rCol", &rCol, def | low,
			     "Half-light radius (in pixels)", 6, 1);
	parameters.addMember("g1Col", &g1Col, def | low,
			     "g1 shape estimate (use native shape if 0)", 7, 0);
	parameters.addMember("g2Col", &g2Col, def | low,
			     "g2 shape estimate (use native shape if 0)", 8, 0);
	parameters.addMember("aCol", &aCol, def | low,
			     "Major axis (in WC)", 7, 1);
	parameters.addMember("bCol", &bCol, def | low,
			     "Minor axis (in WC)", 8, 1);
	parameters.addMember("paCol", &paCol, def | low,
			     "Position angle (in WC)", 9, 1);
	parameters.addMember("fwdColStart", &fwdColStart, def | low,
			     "first forwarded column from input catalog", 0, 1);
	parameters.addMember("fwdColEnd", &fwdColEnd, def | low,
			     "last forwarded column from input catalog", 0, 1);
	parameters.addMember("maxBadPixels", &maxBadPixels, def,
			     "Maximum fraction of bad pixels in postage stamp", 0.1);
	parameters.addMember("maxBadFlux", &maxBadFlux, def,
			     "Maximum fraction of model flux in bad pixels", 0.05);
	parameters.addMember("interpolationOrder", &interpolationOrder, def,
			     "GL order for bad pixel interpolation", 4);
    }
    
    parameters.setDefault();
    
    try {
	for (int i=1; i<argc; i++) {
	    // Open & read all specified input files
	    ifstream ifs(argv[i]);
	    if (!ifs) {
		cerr << "Can't open file " << argv[i] << endl;
		exit(1);
	    }
	    try {
		parameters.setStream(ifs);
	    } catch (std::runtime_error &m) {
		cerr << "In file " << argv[i] << ":" << endl;
		quit(m,1);
	    }
	}
	// ... and from stdin
	try {
	    parameters.setStream(cin);
	} catch (std::runtime_error &m) {
	    cerr << "In stdin:" << endl;
	    quit(m,1);
	}
	
	cerr << "# " << stringstuff::taggedCommandLine(argc, argv) << endl;
	parameters.dump(cerr);
	
	// (1) Read the PSF
	PSFExModel *model;
	try {
	    model = new PSFExModel(psfName.c_str());
	} catch (PSFExError &e) {
	    cerr << "Error reading PSFEx model: " << e.what() << "; exiting." << endl; 
	    exit(1);
	} catch (...) {
	    cerr << "Error reading PSFEx model; exiting." << endl; 
	    exit(1);
	}
	
	if (!model->selfTest(1, 3, psfOrder)) { // this is a very lenient self-test
	    cerr << "PSFEx model did not pass self test; exiting." << endl;
	    exit(1);
	}
	
	cerr << "reading image " << flush;
	// (2) Open image
	Image<> sci;
	FITSImage<> fits(fitsName);
	sci = fits.extract();
        
	cerr << ", header " << flush;
	// (3) Read astrometry
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
	
#ifdef DEBUGFDNTPSFEX
	cerr << "# WCS " << xw << " " << yw << " is the (0,0) pixel WC w.r.t. the TP" << endl;
#endif
	
	cerr << ", weight " << flush;
	// (4) Open weight image and read weight scale
	Image<> wt;
	FITSImage<> wtfits(weightName);
	wt = wtfits.extract();
	cerr << "read, " << flush;
	
	HdrRecordBase* weightScaleRecord = h.find(weightScaleKey);
	double weightScale = 1.0;
	if (weightScaleRecord)
	    weightScale = atof(weightScaleRecord->getValueString().c_str());
	else
	    cerr << "WARNING: weight scale key not found in header; assuming weight scale is 1.0\n";
	
	cerr << "SCAMP weight scale: " << weightScale << endl;
	
	HdrRecordBase* fluxScaleRecord = h.find(fluxScaleKey);
	double fluxScale = 1.0;
	if (fluxScaleRecord)
	    fluxScale = atof(fluxScaleRecord->getValueString().c_str());
	else
	    cerr << "WARNING: flux scale key not found in header; assuming flux scale is 1.0\n";
	
	cerr << "SCAMP flux scale: " << fluxScale << endl;
	
	cerr << "reading segmentation map" << endl;
	Image<> seg;
	FITSImage<> fitsseg(segmentationName.c_str());
	seg = fitsseg.extract();
	
	// (5) rescale weight map
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
	
	// initialize sums
	UnweightedShearEstimator se;
	int ngood=0;
	int nfail=0;
	int nfail_post=0;
	int nfail_postmeas=0;
	int nfail_badpix=0;
	int nfail_bpgl=0;
	int nfail_bpff=0;
	int nfail_fdnt=0;
	int nfail_seg=0;
	
	// (6) Loop through objects
	int nread=MAX(idCol, raCol);
	nread = MAX(nread, decCol);
	nread = MAX(nread, segIdCol);
	nread = MAX(nread, rCol);
	nread = MAX(nread, aCol);
	nread = MAX(nread, bCol);
	nread = MAX(nread, paCol);
	nread = MAX(nread, bgCol);
	nread = MAX(nread, g1Col);
	nread = MAX(nread, g2Col);
	nread = MAX(nread, fwdColEnd);
	vector<string> readvals(nread);
	ifstream ccat(catName.c_str());
	if (!ccat) {
	    cerr << "Error opening catalog file " << catName << endl;
	    exit(1);
	}
	string buffer;
	tmv::SymMatrix<double> covE(2);
	
	cout << "# id x_pix y_pix eta1 eta2 sig1 sig2 cov12 mu egFix fdFlags "
	     << "considered_success forwarded_column" << endl;
	
	
	while (stringstuff::getlineNoComment(ccat, buffer)) { // all objects
	    /*
	    cerr << "######################################################" << endl
		 << "// (a) Acquire info of object " << flush;
	    */
	    istringstream iss(buffer);
	    for (int i=0; i<nread; i++) iss >> readvals[i];
	    if (!iss) {
		cerr << "Bad catalog input line: " << buffer;
		exit(1);
	    }
	    string id = readvals[idCol-1];
	    //cerr << id << endl;
	    double ra0 = atof(readvals[raCol-1].c_str());
	    double dec0 = atof(readvals[decCol-1].c_str());
	    double a_wc = atof(readvals[aCol-1].c_str());
	    double b_wc = atof(readvals[bCol-1].c_str());
	    double pa_wc = atof(readvals[paCol-1].c_str());
	    double r_pix = atof(readvals[rCol-1].c_str()); // FLUX_RADIUS in pixels, not wcs
	    double bg = 0.;
	    if (bgCol > 0)
		bg = atof(readvals[bgCol-1].c_str());
	    string fwd = "";
	    if (fwdColStart <= fwdColEnd && fwdColEnd > 0) {
		for (int i=max(1,fwdColStart); i<=fwdColEnd; i++) {
		    fwd += readvals[i-1];
		    fwd += " ";
		}
	    }
	    
	    double x_wc, y_wc; // tangent plane coordinates
	    double x_pix, y_pix;
	    fullmap->toWorld(x_pix, y_pix, x_wc, y_wc);
	    
#ifdef DEBUGFDNTPSFEX
	    cerr << "(a) processing object at (RA,dec)=(" << ra0 << "," << dec0 << ")=(" 
		 << x_wc << "," << y_wc << ") " << "(x,y)=(" << x_pix << "," << y_pix << ")" << endl;
#endif
	    
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
	    
	    if (!safe_chip_bounds.includes(int(x_pix), int(y_pix)))
		continue;

	    cerr << "######################################################" << endl
		 << "// (a) Acquire info of object " << id << endl
		 << "# object at WCS " << x_wc << " " << y_wc << endl; // DEBUG
#ifdef DEBUGFDNTPSFEX
	    cerr << "processing object at (RA,dec)=(" << ra0 << "," << dec0 << ")=(" 
		 << x_wc << "," << y_wc << ") " << "(x,y)=(" << x_pix << "," << y_pix << ")" << endl;
#endif
	    
#ifdef DEBUGFDNTPSFEX
	    cerr << "// (b) get the Jacobian of coordinate transformation at that position" << endl;
#endif
	    Matrix22 J = fullmap->dWorlddPix(x_pix,y_pix);
	    double rescale = sqrt(fabs(J(0,0)*J(1,1)-J(0,1)*J(1,0)));
	    // distorted approximate pixel scale: how much the scales are enlarged by WC transform
	    
#ifdef DEBUGFDNTPSFEX
	    cerr << "rescaling by " << rescale << endl;
	    
	    cerr << "// (c) get PSF model at the position" << endl;
#endif
	    // generate model PSF at this location
	    model->fieldPosition(x_pix, y_pix); // model now holds the psf model at this position
	    cerr << "flux before normalization: " << model->sb()->getFlux() << endl;
	    model->setFlux(1.0);
	    
	    // FWHM is twice the half-light radius for Gaussian psf, scaled to world coords
	    double ee50psf = model->getFWHM()/2.*rescale; 
	    Ellipse psfBasis(0., 0., log(ee50psf/1.17));
	    
	    // generate PSF image in WC
	    SBDistort psfWCS(*(model->sb()),J(0,0),J(0,1),J(1,0),J(1,1));
	    cerr << "dWorld/dPix " << J;  // J comes with its own endl
	    cerr << ee50psf << endl;
	    // sets flux of basis correctly to have flux normalization after distortion
	    psfWCS.setFlux(1.); 
	    
	    double dx = ee50psf / 2.35; // magic 4.7 factor from PSFEx
	    cerr << "drawing psf postage stamp with dx=" << dx << endl;
	    Image<> ipsf = psfWCS.draw(dx);
	    
	    // Measure PSF GL size & significance
	    Image<> psfwt(ipsf.getBounds());
	    psfwt = pow(1e-5/dx, -2.);
	    psfBasis.setMu( psfBasis.getMu() - log(dx));
	    GLSimple<> gl(ipsf, psfwt, psfBasis, 4);
	    if (!gl.solve()) {
		cerr << "# Failed measuring psf, flags " << gl.getFlags() 
		     << "; ignoring object" << endl;
		continue;
	    }
	    psfBasis = gl.getBasis();
	    psfBasis.setMu( psfBasis.getMu() + log(dx));
#ifdef DEBUGFDNTPSFEX
	    cerr << "# PSF e and GL sigma: " << psfBasis.getS()
		 << " " << exp(psfBasis.getMu())
		 << endl;
#endif
	    
	    // convert shear
	    double g1, g2;
	    psfBasis.getS().getG1G2(g1, g2);

	    cout << id 
		 << " " << ra0
		 << " " << dec0
		 << " " << x_pix
		 << " " << y_pix 
		//<< " " << mag 
		 << " " << g1
		 << " " << g2
		 << " " << sig1
		 << " " << sig2
		 << " " << cov12
		 << " " << psfBasis.getMu()
		 << " " << gl.getFlags()
		 << " " << fwd
		 << endl;
	    
	    delete fep;
	}
	
	delete model;
	delete fullmap;
	
	// Print out mean shear estimate
	Shear S(se);
	double g1, g2, sig1, sig2;
	S.getG1G2(g1,g2);
	se.sigmaE(sig1,sig2, false);
	
	cout << fixed << setprecision(6);
	// Approximate the reduced-shear error as 1/2 of the distortion error
	cout << "# Means: " << g1 << " +- " << sig1/2
	     << " " << g2 << " +- " << sig2/2
	     << endl;
	cout << "# Failed / good " << nfail << " / " << ngood << endl;
	cout << "# Reasons for failure:\n# " << nfail_post 
	     << "\t couldn't get postage stamp\n# " 
	     << nfail_badpix << "\t had too many bad pixels\n# " 
	     << nfail_bpgl << "\t couldn't GL-interpolate bad pixels\n# " 
	     << nfail_bpff << "\t had too much flux in bad pixels\n# " 
	     << nfail_postmeas << "\t failed after measurement\n# "  
	     << nfail_fdnt << "\t had some error in FDNT\n# " 
	     << nfail_seg << "\t had ambiguous segmentation information\n"<< endl;
    } catch (std::runtime_error &m) {
	cerr << m.what() << endl;
	quit(m,1);
    }
}
