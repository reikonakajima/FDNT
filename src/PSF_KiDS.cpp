// PSF_KiDS.cpp
// Run GLSimple on images for the PSFs.  Takes co-add star catalog as input.
//
// input:   [stdin] coadd star catalog, in (RA,Dec) coordinates
// usage:   PSF_KiDS <SCI image> <WT image> <SCAMP pixel map> <segmentation image> <PSF GL Fit Order>
// output:  [stdout] columns of (0) ID (1) RA (2) DEC (3) xpix  (4) ypix  (5) EE50(pix)  (6) mag
//                  (7) g1 (8) g2 (9) sigE1 (10) sigE2 (11) covar (12) mu (13) signif (14) flags
// outline: 1) read in SCI and WT chip FITS images
//          2) read in SCAMP pixel map
//        For each catalog item from stdin:
//          3) find pixel coordinates of each star
//          4) meausure PSF size and shape
//          5) print output
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
    
    // Input files
    string fitsName = argv[1];
    string weightName = argv[2];
    string wcsName = argv[3];
    string segmentationName = argv[4];
    int    segmentationPadding = 0;
    int    order = atoi(argv[5]);
    
    // Data properties
    double sky = 0.;
    int stampSize = 64;
    string weightScaleKey = "WTSCALE";
    
    // GL interpolation of missing values
    double maxBadPixels = 0.1;  // Maximum fraction of bad pixels in postage stamp
    double maxBadFlux = 0.05;  // Maximum fraction of model flux in bad pixels
  
    try {

	cerr << "reading image" << flush;
	// (1) read in SCI and WT chip FITS images
	Image<> sci;
	FITSImage<> fits(fitsName);
	sci = fits.extract();
        
	cerr << ", weight" << flush;
	Image<> wt;
	FITSImage<> wtfits(weightName);
	wt = wtfits.extract();

	// (2) read in astrometry
	cerr << ", header " << flush;
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
	
	//cerr << "# WCS " << xw << " " << yw << " is the (0,0) pixel WC w.r.t. the TP" << endl;
	
	cerr << "read, " << flush;
	
	cerr << "reading segmentation map" << endl;
	Image<> seg;
	FITSImage<> fitsseg(segmentationName.c_str());
	seg = fitsseg.extract();
	
	// get save_chip_bounds
	Bounds<int> chip_bounds = sci.getBounds();
	Bounds<int> safe_chip_bounds = chip_bounds;
	safe_chip_bounds.addBorder(-2);
	
	// initialize sums
	int ngood = 0;
	int nfail = 0;
	int nfail_seg = 0;
	int nfail_post = 0;
	int nfail_badpix = 0;
	int nfail_bpgl = 0;
	int nfail_bpff = 0;
	int nfail_sig = 0;
	
        //  Coadd (star/galaxy) catalog columns [from stdin]
	//
        //  1   1 FIELD_POS       Reference number to field parameters
        //  2   2 SeqNr           Running object number
        //  3 113 MAG_BEST        Best of MAG_AUTO and MAG_ISOCOR                 [mag]
        //  4 116 BackGr          Background at centroid position                 [count]
        //  5 117 Level           Detection threshold above background            [count]
        //  6 133 ALPHA_J2000     Right ascension of barycenter (J2000)           [deg]
        //  7 134 DELTA_J2000     Declination of barycenter (J2000)               [deg]
        //  8 138 X2_WORLD        Variance along X-WORLD (alpha)                  [deg**2]
        //  9 139 Y2_WORLD        Variance along Y-WORLD (delta)                  [deg**2]
        // 10 140 XY_WORLD        Covariance between X-WORLD and Y-WORLD          [deg**2]
        // 11 149 A_WORLD         Profile RMS along major axis (world units)      [deg]
        // 12 150 B_WORLD         Profile RMS along minor axis (world units)      [deg]
        // 13 152 THETA_WORLD     Position angle (CCW/world-x)                    [deg]
        // 14 156 ELLIPTICITY     1 - B_IMAGE/A_IMAGE
        // 15 178 FWHM_WORLD      FWHM assuming a gaussian core                   [deg]
        // 16 187 Flag            Extraction flags
        // 17 188 FLUX_RADIUS     Fraction-of-light radii                         [pixel]
        // 18 193 MASK            mask value (addmask)            0==not in masked region

	// define catalog item columns (0-indexed)
	int idCol=2, raCol=6, decCol=7, aCol=11, bCol=12, paCol=13, rCol=17, bgCol=4, magCol=3;
	int nread=18;
	vector<string> readvals(nread);

	// for each catalog item
	string buffer;
	cout << "# id   ra   dec  xpix  ypix  EE50(pix)  mag  g1   g2  mu  "
	     << "signif  model_flux  flags" << endl;
	while (stringstuff::getlineNoComment(cin, buffer)) {

	    /*/ // DEBUG messages
	    cerr << "######################################################" << endl
		 << "// (a) Acquire info of object " << flush;
	    /**/
	    istringstream iss(buffer);
	    for (int i=0; i<nread; i++) iss >> readvals[i];
	    if (!iss) {
		cerr << "Bad catalog input line: " << buffer;
		exit(1);
	    }
	    string id = readvals[idCol-1];
	    //cerr << "ID: " << id << endl; // DEBUG
	    double ra0 = atof(readvals[raCol-1].c_str());
	    double dec0 = atof(readvals[decCol-1].c_str());
	    double a_wc = atof(readvals[aCol-1].c_str());
	    double b_wc = atof(readvals[bCol-1].c_str());
	    double pa_wc = atof(readvals[paCol-1].c_str());
	    double r_pix = atof(readvals[rCol-1].c_str()); // FLUX_RADIUS in pixels, not wcs
	    double bg = atof(readvals[bgCol-1].c_str());
	    double mag = atof(readvals[magCol-1].c_str());

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

	    /*/ // DEBUG MESSAGES
	    cerr << "######################################################" << endl
		 << "// (a) Acquire info of object " << id << endl
		 << "# object at WCS " << x_wc << " " << y_wc << endl; // DEBUG
	    cerr << "processing object at (RA,dec)=(" << ra0 << "," << dec0 << ")=(" 
		 << x_wc << "," << y_wc << ") " << "(x,y)=(" << x_pix << "," << y_pix << ")" << endl;
	    /**/
	    
	    //cerr << "// (b) Starting ellipse of star" << endl;   // DEBUG

	    double e1start = (a_wc*a_wc-b_wc*b_wc)/(a_wc*a_wc+b_wc*b_wc);
	    double e2start = e1start * sin(2*pa_wc*PI/180.);
	    e1start *= cos(2*pa_wc*PI/180.);
	    
	    //cerr << "// (c) building FitExposure" << endl;  // DEBUG
	    FitExposure<>* fep;

	    //cerr << "// get the postage stamp" << endl;  // DEBUG
	    double xpix_stamp = 0.;  // postage stamp centers
	    double ypix_stamp = 0.;

	    int xi0 = x_pix-stampSize/2;
	    if (xi0 < chip_bounds.getXMin()) {
		xpix_stamp = xi0 - 1;
		xi0 = chip_bounds.getXMin();
	    }
	    if ((xi0 + stampSize+1) > chip_bounds.getXMax()) 
		xi0 = chip_bounds.getXMax() - stampSize - 1;
	    int yi0 = y_pix-stampSize / 2;
	    if (yi0 < chip_bounds.getYMin()) {
		ypix_stamp = yi0 - 1;
		yi0 = chip_bounds.getYMin();
	    }
	    if ((yi0 + stampSize+1) > chip_bounds.getYMax()) 
		yi0 = chip_bounds.getYMax() - stampSize - 1;
	    /*/ // DEBUG
	    cerr << "// " << xi0 << " " << xi0+stampSize-1 << " " 
		 << yi0 << " " << yi0+stampSize-1 << endl;
	    cerr << "// hopefully inside x=" << chip_bounds.getXMin() 
		 << ".." << chip_bounds.getXMax() 
		 << ", y=" << chip_bounds.getYMin() 
		 << ".." << chip_bounds.getYMax() << endl; 
	    /**/
	    Bounds<int> stamp(xi0, xi0+stampSize-1, yi0, yi0+stampSize-1);
	    //cerr << "set bounds" << endl;
	    Image<float> scistamp = sci.subimage(stamp).duplicate();
	    //cerr << "sci bounds" << endl;
	    Image<float> wtstamp  = wt.subimage(stamp).duplicate();
	    //cerr << "wt  bounds" << endl;
	    Image<float> segstamp = seg.subimage(stamp).duplicate();
	    //cerr << "seg bounds" << endl;

	    xpix_stamp += (scistamp.getBounds().getXMax() + scistamp.getBounds().getXMin()) / 2.;
	    ypix_stamp += (scistamp.getBounds().getYMax() + scistamp.getBounds().getYMin()) / 2.;

	    // all in pure pixel coordinates
	    Ellipse sexE(e1start, e2start, log(r_pix), xpix_stamp, ypix_stamp); 
            /*/  // DEBUG
	    cerr << "starting ellipse r = " << r_pix
		 << "(pix), e = " << e1start << ", " << e2start << endl;
	    /**/	    

	    bool badSegID=false;
	    int segId = 0;
	    // check a small square region
	    for (int dx=0; dx<=1; dx++) {
		for (int dy=0; dy<=1; dy++) {
		    //cerr << int(x_pix) + dx << " " << int(y_pix) + dy << endl;  // DEBUG
		    int dsegId = seg(int(x_pix)+dx, int(y_pix)+dy);
		    if (!segId && dsegId) {  // for the first occurence
			segId = dsegId;
		    }
		    else if (dsegId && segId!=dsegId) {
			badSegID = true;
		    }
		}
	    }
		
	    if (badSegID) {
		cerr << "# fail: could not get segmentation id for " << id << endl;
		/*/ DEBUG: CHECKPLOTS
		  ipsf.shift(1,1);
		  FITSImage<>::writeToFITS("check_"+id+"_psf_badseg.fits",ipsf);
		  wtstamp.shift(1,1);
		  FITSImage<>::writeToFITS("check_"+id+"_wt_badseg.fits",wtstamp);
		  scistamp.shift(1,1);
		  FITSImage<>::writeToFITS("check_"+id+"_sci_badseg.fits",scistamp);
		  segstamp.shift(1,1);
		  FITSImage<>::writeToFITS("check_"+id+"_seg_badseg.fits",segstamp);
		  /*/
		nfail++;
		nfail_seg++;
		continue; 
	    }
	    //cerr << "got segid" << endl;
	    
	    try {
		for (int iy=stamp.getYMin(); iy<=stamp.getYMax(); iy++)
		{
		    for (int ix=stamp.getXMin(); ix<=stamp.getXMax(); ix++) {
			double dxw,dyw;
			fullmap->toWorld(ix,iy,dxw,dyw);
			if (segstamp(ix,iy) && segstamp(ix,iy)!=segId) { //it's another object!
			    for (int iiy=max(stamp.getYMin(),iy-segmentationPadding); 
				 iiy<=min(stamp.getYMax(),iy+segmentationPadding); iiy++) 
				for (int iix=max(stamp.getXMin(),ix-segmentationPadding); 
				     iix<=min(stamp.getXMax(),ix+segmentationPadding); iix++) 
				    wtstamp(iix,iiy)=0.;
			}
		    }
		}
		fep = new FitExposure<>(scistamp, wtstamp, 0, sky+bg);  // in pixel coordinates!
	    }
	    catch (...)
	    {
		cerr << "# fail: could not get postage stamp for " << id << endl;
		nfail++;
		nfail_post++;
		delete fep;
		continue;
	    }

	    //cerr << "// (f) checking bad pixels" << endl;
	    double meanweight=0.;
	    vector< Position<int> > bp = fep->badPixels(meanweight);
	    //cerr << "mean weight: " << meanweight << endl;  // DEBUG
	    if (bp.size() > maxBadPixels*stampSize*stampSize) {
		cerr << "# fail: " << id << " has too many bad pixels" << endl;
		/*/ DEBUG: CHECKPLOTS
		ipsf.shift(1,1);
		FITSImage<>::writeToFITS("check_"+id+"_psf_badpix.fits",ipsf);
		wtstamp.shift(1,1);
		FITSImage<>::writeToFITS("check_"+id+"_wt_badpix.fits",wtstamp);
		scistamp.shift(1,1);
		FITSImage<>::writeToFITS("check_"+id+"_sci_badpix.fits",scistamp);
		segstamp.shift(1,1);
		FITSImage<>::writeToFITS("check_"+id+"_seg_badpix.fits",segstamp);
                /*/
		nfail++;
		nfail_badpix++;
		delete fep;
		continue;
	    }

	    /*/ residual image (science image with the GL model subtracted)
	    Image<float> glstamp = sci.subimage(stamp).duplicate(); 
	    /*/
	    
	    // do GL fit
	    /*/ // DEBUG
	    cerr << "# bad pixel count = " << bp.size() << endl;
	    cerr << "# start ellipse" << sexE << " " <<  order << endl;
	    /**/
	    GLSimple<> star(*fep, sexE, order); 
	    if (!star.solve()) {
		cerr << "# fail: GL fit failed for " << id << " : "
		     << "fail flag: " << star.getFlags() << endl;
		/*/ //DEBUG: CHECKPLOTS
		wtstamp.shift(1,1);
		FITSImage<>::writeToFITS("check_"+id+"_wt_bpgl.fits",wtstamp);
		scistamp.shift(1,1);
		FITSImage<>::writeToFITS("check_"+id+"_sci_bpgl.fits",scistamp);
		segstamp.shift(1,1);
		FITSImage<>::writeToFITS("check_"+id+"_seg_bpgl.fits",segstamp);
		/**/
		nfail++;
		nfail_bpgl++;
		delete fep;
		
		continue;
	    }
	    LVector bvec = star.getB();
	    double flux = bvec.flux();
	    Ellipse basis = star.getBasis();
	    double missingFlux=0.;
	    double scaleFactor = exp(basis.getMu());
	    
	    for (vector< Position<int> >::iterator it=bp.begin(); it<bp.end(); it++) {
		LVector psi(bvec.getOrder());
		Position<double> xunit = basis.inv(Position<double>((*it).x, (*it).y));
		psi.fillBasis(xunit.x, xunit.y, scaleFactor);
		scistamp((*it).x,(*it).y)=bvec.dot(psi);
		wtstamp((*it).x,(*it).y)=meanweight;
		// The interpolated pixel is assumed to have the mean weight of the good pixels;
		// this makes sense because in the Fourier code homogeneous uncertainties are 
		// assumed
		missingFlux += scistamp((*it).x,(*it).y);
		scistamp((*it).x,(*it).y)+=sky+bg;
	    }
	    //missingFlux *= rescale*rescale; // scale flux with coordinate system
	    
	    /*/ DEBUG
	      for (int iy=stamp.getYMin(); iy<=stamp.getYMax(); iy++) {
	      for (int ix=stamp.getXMin(); ix<=stamp.getXMax(); ix++) {
	      LVector psi(bvec.getOrder());
	      Position<double> xunit = 
	      basis.inv(Position<double>(xwstamp(ix,iy), ywstamp(ix,iy)));
	      psi.fillBasis(xunit.x, xunit.y, scaleFactor);
	      glstamp(ix,iy) -= bvec.dot(psi);   
	      }
	      }
	      
	      cerr << "# flux = " << flux << endl;
	      cerr << "# scaleFactor = " << scaleFactor << endl;
	      cerr << "# missingFlux = " << missingFlux << endl;
	      /*/
	    
	    if (missingFlux/flux > maxBadFlux) {
		cerr << "# fail: bad pixels in " << id << " have too high flux fraction" << endl;
		/*/ DEBUG: CHECKPLOTS
		  ipsf.shift(1,1);
		  FITSImage<>::writeToFITS("check_"+id+"_psf_bpff.fits",ipsf);
		  wtstamp.shift(1,1);
		  FITSImage<>::writeToFITS("check_"+id+"_wt_bpff.fits",wtstamp);
		  scistamp.shift(1,1);
		  FITSImage<>::writeToFITS("check_"+id+"_sci_bpff.fits",scistamp);
		  segstamp.shift(1,1);
		  FITSImage<>::writeToFITS("check_"+id+"_seg_bpff.fits",segstamp);
		  glstamp.shift(1,1);
		  FITSImage<>::writeToFITS("check_"+id+"_glr_bpff.fits",glstamp);
		  /*/
		nfail++;
		nfail_bpff++;
		delete fep;
		continue;
	    }
	    
	    // get shape info
	    double g1, g2;
	    basis.getS().getG1G2(g1, g2);
	    /*/
	    tmv::SymMatrix<double> cov;  // unavailable
	    cov = gl.covar();
	    double varE1, varE2, covar;
	    varE1 = cov(0,0);
	    varE2 = cov(1,1);
	    covar = cov(1,0);
            /*/
	    double significance, b00, var_b00;
	    star.b00(b00, var_b00);
	    if (b00 == -1)
		significance = -1.;
	    else
		significance = b00 / sqrt(var_b00);

	    if (significance < 100) {
		cerr << "# fail: "<< id << " has S/N < 100" << endl;
		nfail++;
		nfail_sig++;
		continue;
	    }

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
		 << fixed << setprecision(5) << setw(8)
		 << g1 << " " << setw(8)
		 << -g2 << " "  // flip g2 so that it matches ra/dec (not xpix/ypix)
		 << basis.getMu() << " ";
	    cout.flags(old_settings);
	    cout << setw(9)
		 << significance << " " << setw(6)
		 << flux << " "
		 << star.getFlags()
		 << endl;

	    ++ngood;
	    delete fep;
	}
	
	delete fullmap;
	
	cerr << "# Failed / good " << nfail << " / " << ngood << endl;
	cout << "# Failed / good " << nfail << " / " << ngood << endl;
	cout << "# Reasons for failure:\n# " << nfail_post 
	     << "\t couldn't get postage stamp\n# " 
	     << nfail_badpix << "\t had too many bad pixels\n# " 
	     << nfail_bpgl << "\t GL-fit failed\n# " 
	     << nfail_bpff << "\t had too much flux in bad pixels\n# " 
	     << nfail_sig << "\t failed the significance test\n# " 
	     << nfail_seg << "\t had ambiguous segmentation information\n";
    } catch (std::runtime_error &m) {
	cerr << m.what() << endl;
	quit(m,1);
    }
}
