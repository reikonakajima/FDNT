// PSFEx_KiDS.cpp
// Generate PSFEx psf models and measure shape; generate catalog
//
// input:   [stdin] coadd star catalog, in (RA,Dec) coordinates
// usage:   PSFEx_KiDS <PSFEx file> <pixel map (SCAMP header files)> <PSF GL Fit Order>
// output:  [stdout] columns of (0) ID (1) RA (2) DEC (3) xpix  (4) ypix  (5) EE50(pix)  (6) mag
//                   (7) g1 (8) g2 (9) sigE1 (10) sigE2 (11) covar (12) mu (13) signif (14) flags
// outline: 1) read in PSF file
//          2) read in pixel map file
//        For each catalog item from stdin:
//          3) find pixel coordinates of each star
//          4) get PSF model at that pixel coordinate
//          5) meausure PSF model size and shape
//          6) print output
// 

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
    
    // The fit:
    double maskSigma;
    int    psfOrder = atoi(argv[4]);
    
    // Input files
    string psfName = argv[1];
    int chipNumber = atoi(argv[2]);
    string wcsName = argv[3];
    
    // Data properties
    double sky = 0.;
    //double rn = 1000.;
    //double gain = 1.;
    int stampSize = 64;
    
    try {

	// (1) Read the PSF
	PSFExModel *model;
	try {
	    model = new PSFExModel(psfName.c_str(), chipNumber);
	} catch (PSFExError &e) {
	    cerr << "Error reading PSFEx model: " << e.what() << "; exiting." << endl; 
	    exit(1);
	}
	
	if (!model->selfTest(1, 3, psfOrder)) { // this is a very lenient self-test
	    cerr << "PSFEx model did not pass self test; exiting." << endl;
	    exit(1);
	}
	
	// (2) Read astrometry
	img::ImageHeader h;
	ifstream mapfs(wcsName.c_str());
	Image<> sci;
	if (!mapfs) {
	    cerr << "Could not open astrometry file " << wcsName
		 << "; trying FITS header in science frame" << endl;
	    try {
		FITSImage<> fits(wcsName);
		sci = fits.extract();
	    } catch (AstrometryError& e) {
		cerr << "Could not get header info from science (coadd) frame" << endl;
		cerr << e.what() << endl;
		exit(0);
	    }
	    h = *(sci.header());
	}
	else {
	    h = img::HeaderFromStream(mapfs);
	}
	SCAMPMap *fullmap = new SCAMPMap(h,0);
	double xw, yw;
	fullmap->toWorld(0,0,xw,yw);

	// define chip bounds
        Bounds<int> chip_bounds(1,2040,1,4050);  // I hope it's the same for every chip...
        Bounds<int> safe_chip_bounds = chip_bounds;
        safe_chip_bounds.addBorder(-2);

	// track fits
	int ngood = 0;
	int nfail = 0;

	// define catalog item columns (0-indexed)
	int raCol=0, decCol=1, magCol=2;
	int nread=3;
	vector<string> readvals(nread);

	// For each catalog item
	string buffer;
	cout << "# id   ra   dec  xpix  ypix  EE50(pix)  mag  g1   g2  mu  "
	     << "signif  model_flux  flags" << endl;
	while (stringstuff::getlineNoComment(cin, buffer)) {

	    /*/ DEBUG messages
	    cerr << "######################################################" << endl
		 << "// (a) Acquire info of object " << flush;
	    /*/
	    istringstream iss(buffer);
	    for (int i=0; i<nread; i++) iss >> readvals[i];
	    if (!iss) {
		cerr << "Bad catalog input line: " << buffer;
		exit(1);
	    }
	    string id = "00";
	    // cerr << "ID: " << id << endl; // DEBUG
	    double ra0 = atof(readvals[raCol].c_str());
	    double dec0 = atof(readvals[decCol].c_str());
	    double r_pix = -9.9; // FLUX_RADIUS in pixels, not wcs
	    double mag = atof(readvals[magCol].c_str());
	    
	    double x_wc, y_wc; // tangent plane coordinates
	    double x_pix, y_pix;
	    
	    // deproject spherical coordinates (ra0, dec0) to xi (x_wc) and eta (y_wc)
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
		cerr << "(b) processing object at (RA,dec)=(" << ra0 << "," << dec0 << ") => " 
		     << "(xi,eta)=(" << x_wc << "," << y_wc << ")" << endl;
		cerr << "toPix failure" << endl;
		exit(0);
	    }

	    if (!safe_chip_bounds.includes(int(x_pix), int(y_pix)))
		continue;

            // (b) get the Jacobian of coordinate transformation at that position
	    Matrix22 J = fullmap->dWorlddPix(x_pix,y_pix);
	    double rescale = sqrt(fabs(J(0,0)*J(1,1)-J(0,1)*J(1,0)));
	    // distorted approximate pixel scale: how much the scales are enlarged by WC transform
	    
	    // (c) get PSF model at the position
	    model->fieldPosition(x_pix, y_pix); // model now holds the psf model at this position
	    double original_flux = model->sb()->getFlux(); // observe how the model flux changes
	    model->setFlux(1.0);
	    
	    // FWHM is twice the half-light radius for Gaussian psf, scaled to world coords
	    double ee50psf = model->getFWHM()/2.*rescale; 
	    Ellipse psfBasis(0., 0., log(ee50psf/1.17));
	    
	    // generate PSF image in WC
	    SBDistort psfWCS(*(model->sb()),J(0,0),J(0,1),J(1,0),J(1,1));
	    // sets flux of basis correctly to have flux normalization after distortion
	    psfWCS.setFlux(1.); 
	    
	    double dx = ee50psf / 2.35; // psfWCS: magic 4.7 factor from PSFEx
	    //cerr << "drawing psf postage stamp with dx=" << dx << endl;
	    Image<> ipsf = psfWCS.draw(dx);
	    
	    // Measure PSF GL size & significance
	    Image<> psfwt(ipsf.getBounds());
	    psfwt = pow(1e-5/dx, -2.);
	    psfBasis.setMu( psfBasis.getMu() - log(dx));
	    GLSimple<> gl(ipsf, psfwt, psfBasis, 4);
	    if (!gl.solve()) {
		cerr << "# Failed measuring psf, flags " << gl.getFlags() 
		     << "; ignoring object " << ra0 << " " << dec0 << endl;
		nfail += 1;
		continue;
	    }
	    psfBasis = gl.getBasis();
	    psfBasis.setMu( psfBasis.getMu() + log(dx));
	    ngood += 1;

	    /*/ DEBUG message
	    cerr << "# PSF e and GL sigma: " << psfBasis.getS()
		 << " " << exp(psfBasis.getMu())
		 << endl;
	    /*/
	    /*/ DEBUG plot
	    ipsf.shift(1,1);
	    FITSImage<>::writeToFITS("check_"+id+"_psfex.fits",ipsf);
	    /*/
	    // get shape info
	    double g1, g2;
	    psfBasis.getS().getG1G2(g1, g2);
	    /*/
	    tmv::SymMatrix<double> cov;  // unavailable
	    cov = gl.covar();
	    double varE1, varE2, covar;
	    varE1 = cov(0,0);
	    varE2 = cov(1,1);
	    covar = cov(1,0);
            /*/
	    double significance, b00, var_b00;
	    gl.b00(b00, var_b00);
	    if (b00 == -1)
		significance = -1.;
	    else
		significance = b00 / sqrt(var_b00);

	    ios::fmtflags old_settings = cout.flags();
	    int old_precision = cout.precision();
	    cout << fixed << setprecision(6);
	    cout << id << " "
		 << fixed << setprecision(6)
		 << ra0 << " "
		 << dec0 << " "
		 << fixed << setprecision(2) << setw(7)
		 << x_pix << " " << setw(7)
		 << y_pix << " "
		 << fixed << setprecision(3) << setw(5)
		 << r_pix << " "
		 << fixed << setprecision(2) << setw(5)
		 << mag << " "
		 << fixed << setprecision(5)
		 << g1 << " "
		 << g2 << " "
		 << psfBasis.getMu() << " ";
	    cout.flags(old_settings);
	    cout << setprecision(4) << setw(9)
		 << significance << " " << setw(5)
		 << original_flux << " "
		 << gl.getFlags()
		 << endl;
	}
	
	delete model;
	delete fullmap;
	cout << "# Failed / good " << nfail << " / " << ngood << endl;

    } catch (std::runtime_error &m) {
	cerr << m.what() << endl;
	quit(m,1);
    }
}
