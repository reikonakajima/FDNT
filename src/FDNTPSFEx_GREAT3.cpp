// FDNTPSFEx_GREAT3.cpp
// Run FDNT on images with PSFEx psf model and SExtractor input catalog

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
#include "PolyMap.h"

#define DEBUGFDNTPSFEX
//#define CHECKPLOTS

using namespace laguerre;
using namespace sbp;
using namespace astrometry;

int
main(int argc, 
     char *argv[])
{
  // Parameters to set:
  // columns for necessary input:
  int idCol;
  int segIdCol;
  int raCol;
  int decCol;
  //int magCol;
  int bgCol;
  int rCol;
  int aCol,bCol,paCol; // WC observed ellipse
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
      parameters.addMember("fitsName",&fitsName, def,
			   "Image FITS file", "image-000-0.fits");
      parameters.addMember("catName",&catName, def,
			   "Input SExtractor catalog", "sextractor-000-0.cat");
      parameters.addMember("psfName",&psfName, def,
			   "Input PSFEx file", "psfex-000-0.psf");
      parameters.addMember("psfOrder",&psfOrder, def | low,
			   "Maximum order of polynomial psf variation used; -1 for order of PSFEx model", -1, -1);
      parameters.addMember("weightName",&weightName, def, "Input weight map proportional to inverse variance", "weight-000-0.fits");
      parameters.addMember("wcsName",&wcsName, def,
			   "World coordinate header file", "wcs_is_pixel.txt");
      parameters.addMember("weightScaleKey",&weightScaleKey, def,
			   "Scaling factor for weight map (e.g. from Swarp) keyword in WCS header", "WTSCALE");
			   // this is the scale that transforms weight proportional to 1/sig^2 to weights EQUAL to 1/sig^2 of the respective single frame
			   // note that single frame rescaling is happening on top of this on both the science and weight frames
      parameters.addMember("fluxScaleKey",&fluxScaleKey, def,
			   "Scaling factor for flux (e.g. from WCS header) keyword in WCS header", "FLXSCALE");
      parameters.addMember("segmentationName",&segmentationName, def,
			   "segmentation map", "segmentation-000-0.fits");
      parameters.addMember("segmentationPadding",&segmentationPadding, def | low,
			   "padding around segmentation maps [pix]", 0, 0);
      parameters.addMember("stampSize",&stampSize, def | low,
			   "Pixel size of img postage stamps", 48, 3);
      parameters.addMemberNoValue("DATA PROPERTIES:",0,
				  "Input data");
      parameters.addMember("sky",&sky, def,
			   "sky level", 0.);
      //parameters.addMember("rn",&rn, def | low,
	//		   "Read noise, i.e. RMS at 0 ADU", 1000., 0.);
      //parameters.addMember("gain",&gain, def | low,
	//		   "Gain, e/ADU for Poissson, 0 for no Poisson", 1., 0.);
      parameters.addMemberNoValue("FDNT:",0,
				  "Fitting parameters:");
      parameters.addMember("order",&order, def,
			   "FDNT GL-fit order (0 for FFT)", 0);
      parameters.addMember("maskSigma",&maskSigma, def,
			   "GL and FDNT mask sigma (-1 auto)", -1. /*4.*/); // changed 21.11.12
      parameters.addMemberNoValue("CATALOG TRANSLATION:",0,
				  "Columns holding input data");
      parameters.addMember("idCol",&idCol, def | low,
			   "Object ID", 1, 1);
      parameters.addMember("segIdCol",&segIdCol, def | low,
			   "Object ID in segmentation map; determined automatically if 0", 0, 0);
      parameters.addMember("raCol",&raCol, def | low,
			   "RA centroid", 2, 1);
      parameters.addMember("decCol",&decCol, def | low,
			   "DEC centroid", 3, 1);
      //parameters.addMember("magCol",&magCol, def | low,
	//		   "Magnitude", 4, 1);
      parameters.addMember("bgCol",&bgCol, def | low,
			   "photometric background flux (zero if bgCol==0)", 0, 0);
      parameters.addMember("rCol",&rCol, def | low,
			   "Half-light radius (in pixels)", 6, 1);
      parameters.addMember("g1Col",&g1Col, def | low,
			   "g1 shape estimate (pixel coordinates)", 0, 0);  // use native (g1,g2)
      parameters.addMember("g2Col",&g2Col, def | low,
			   "g2 shape estimate (pixel coordinates)", 0, 0);
      parameters.addMember("aCol",&aCol, def | low,
			   "Major axis (in WC)", 7, 1);
      parameters.addMember("bCol",&bCol, def | low,
			   "Minor axis (in WC)", 8, 1);
      parameters.addMember("paCol",&paCol, def | low,
			   "Position angle (in WC)", 9, 1);
      parameters.addMember("fwdColStart",&fwdColStart, def | low,
			   "first forwarded column from input catalog", 0, 0);
      parameters.addMember("fwdColEnd",&fwdColEnd, def | low,
			   "last forwarded column from input catalog", 0, 0);
      parameters.addMember("maxBadPixels",&maxBadPixels, def,
			   "Maximum fraction of bad pixels in postage stamp", 0.1);
      parameters.addMember("maxBadFlux",&maxBadFlux, def,
			   "Maximum fraction of model flux in bad pixels", 0.05);
      parameters.addMember("interpolationOrder",&interpolationOrder, def,
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
    // and stdin:
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
    try
    {
	model = new PSFExModel(psfName.c_str());
    }
    catch (PSFExError &e)
    {
      cerr << "Error reading PSFEx model: " << e.what() << "; exiting." << endl; 
      exit(1);
    }
    catch (...)
    {
      cerr << "Error reading PSFEx model; exiting." << endl; 
      exit(1);
    }
    
    if (!model->selfTest(1,3,psfOrder)) // this is a very lenient self-test
    {
      cerr << "PSFEx model did not pass self test; exiting." << endl;
      exit(1);
    }
    
    cerr << "reading image " << flush;
    // (3) Open image
    Image<> sci;
    FITSImage<> fits(fitsName);
    sci = fits.extract();
        
    // (2) Read astrometry: use pixel scales!  LinearMap
    img::ImageHeader h;
    cerr << ", header " << flush;
    ifstream mapfs(wcsName.c_str());
    if (!mapfs) {
	cerr << "FIRST" << endl;
	cerr << "Could not open map spec file " << wcsName 
	     << "; trying FITS header in science frame" << endl;
	h = *(sci.header());
    }
    else {
	cerr << "SECOND" << endl;
        h = img::HeaderFromStream(mapfs);
    }
    // setup identity map
    DVector v(6);
    v[0] = 0.; v[1] = 1.; v[2] = 0.; v[3] = 0; v[4] = 0; v[5] = 1.;
    LinearMap *fullmap = new LinearMap(v);
    /*
    double xw,yw;
    fullmap->toWorld(0,0,xw,yw);
#ifdef DEBUGFDNTPSFEX
    cerr << "# WCS " << xw << " " << yw << " is the (0,0) pixel WC w.r.t. the TP" << endl;
#endif
    */
    cerr << ", weight " << flush;

    // (4) Open weight image and read weight scale
    Image<> wt;
    FITSImage<> wtfits(weightName);
    wt = wtfits.extract();
    
    cerr << "read, " << flush;
    
    HdrRecordBase* weightScaleRecord = h.find(weightScaleKey);
    double weightScale = 1.0;
    if (weightScaleRecord) {
	weightScale = atof(weightScaleRecord->getValueString().c_str());
    }
    else
    {
      cerr << "WARNING: weight scale key not found in header; assuming weight scale is 1.0" << endl;
    }
    cerr << "weight scale: " << weightScale << endl;
    
    HdrRecordBase* fluxScaleRecord = h.find(fluxScaleKey);
    double fluxScale = 1.0;
    if (fluxScaleRecord)
	fluxScale = atof(fluxScaleRecord->getValueString().c_str());
    else
	cerr << "WARNING: flux scale key not found in header; assuming flux scale is 1.0" << endl;
    
    cerr << "flux scale: " << fluxScale << endl;
    
    cerr << " and segmentation map" << endl;
    
    // cerr << "reading segmentation map" << endl;
    Image<> seg;
    FITSImage<> fitsseg(segmentationName.c_str());
    seg = fitsseg.extract();
      
    // (5) rescale weight map
    Bounds<int> bounds=sci.getBounds();
    
    weightScale = weightScale/(fluxScale*fluxScale); // that's how the inverse variance scales with flux
    for (int iy=bounds.getYMin(); iy<=bounds.getYMax(); iy++)
	for (int ix=bounds.getXMin(); ix<=bounds.getXMax(); ix++) {
	
	    // weight map just rescaled sky inverse-variance:
	    wt(ix,iy)  = wt(ix,iy)*weightScale;
	    sci(ix,iy) = sci(ix,iy)*fluxScale;
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

    cerr << "NREAD: " << nread << endl;

    vector<string> readvals(nread);
    ifstream ccat(catName.c_str());
    if (!ccat) {
	cerr << "Error opening catalog file " << catName << endl;
	exit(1);
    }
    string buffer;
    tmv::SymMatrix<double> covE(2);

    cout << "# id x_pix y_pix eta1 eta2 sig1 sig2 mu egFix fdFlags considered_success forwarded_column" << endl;


    while (stringstuff::getlineNoComment(ccat, buffer)) { // all objects

	//
	// acquire info of object
	//
      cerr << "######################################################" << endl
	   << "// (a) Acquire info of object " << flush;
      istringstream iss(buffer);
      for (int i=0; i<nread; i++) iss >> readvals[i];
      if (!iss) {
	  cerr << "Bad catalog input line: " << buffer;
	  exit(1);
      }
      string id = readvals[idCol-1];
      cerr << id << endl;
      double ra0 = atof(readvals[raCol-1].c_str());
      double dec0 = atof(readvals[decCol-1].c_str());
      double a_wc = atof(readvals[aCol-1].c_str());
      double b_wc = atof(readvals[bCol-1].c_str());
      double pa_wc = atof(readvals[paCol-1].c_str());
      // FLUX_RADIUS is in pixels! there is no _WORLD alternative
      double r_pix = atof(readvals[rCol-1].c_str());

      // set-up initial g1, g2 from native SExtractor value
      double e1start = (a_wc*a_wc-b_wc*b_wc)/(a_wc*a_wc+b_wc*b_wc);
      double e2start = e1start * sin(2*pa_wc*PI/180.);
      e1start *= cos(2*pa_wc*PI/180.);
      Shear  sexS(e1start, e2start);
      double g1_wc, g2_wc;
      sexS.getG1G2(g1_wc, g2_wc);
      if (g1Col > 0)
	  g1_wc = atof(readvals[g1Col-1].c_str());
      if (g2Col > 0)
	  g2_wc = atof(readvals[g2Col-1].c_str());
      double gabs = sqrt(g1_wc*g1_wc+g2_wc*g2_wc);
      if (gabs > 0.95) { // sanitize shear
	  g1_wc=g1_wc*0.95/gabs;
	  g2_wc=g2_wc*0.95/gabs;
      }

      double bg = 0.;
      if (bgCol > 0)
	  bg = atof(readvals[bgCol-1].c_str());
      string fwd = "";
      if (fwdColStart <= fwdColEnd && fwdColEnd > 0)
      {
	  for (int i=max(1,fwdColStart); i<=fwdColEnd; i++)
	  {
	      fwd += readvals[i-1];
	      fwd += " ";
	  }
      }

      double x_wc = ra0;
      double y_wc = dec0;
      double x_pix, y_pix;
      
      fullmap->toPix(x_wc,y_wc,x_pix,y_pix);
      
#ifdef DEBUGFDNTPSFEX
      cerr << "processing object at (RA,dec)=(" << ra0 << "," << dec0 << ")=(" << x_wc << "," << y_wc << ") " << "(x,y)=(" << x_pix << "," << y_pix << ")" << endl;
#endif
     
#ifdef DEBUGFDNTPSFEX
      cerr << "// (b) get the Jacobian of coordinate transformation at that position" << endl;
#endif
      Matrix22 J = fullmap->dWorlddPix(x_pix,y_pix);
      double rescale = sqrt(fabs(J(0,0)*J(1,1)-J(0,1)*J(1,0)));
      // distorted approximate pixel scale: how much the scales are enlarged by
      // world coordinate transformation
#ifdef DEBUGFDNTPSFEX
      cerr << "rescaling by " << rescale << endl;
#endif

      //
      // get PSF model at the position
      //
#ifdef DEBUGFDNTPSFEX
      cerr << "// (c) get PSF model at the position" << endl;
#endif
      model->fieldPosition(x_pix,y_pix); // model now holds the psf model at this position
      cerr << "flux before normalization: " << model->sb()->getFlux() << endl;
      model->setFlux(1.0);
      // FWHM is twice the half-light radius for Gaussian psf, scaled to world coords
      double ee50psf = model->getFWHM()/2.*rescale;
      Ellipse psfBasis(0., 0., log(ee50psf/1.17));
      SBDistort psfWCS(*(model->sb()),J(0,0),J(0,1),J(1,0),J(1,1));
      cerr << "dWorld/dPix " << J << endl; 
      cerr << ee50psf << endl;
      // set flux of basis correctly to have flux normalization after distortion
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
	  cerr << "# Failed measuring psf, flags " << gl.getFlags() << "; ignoring object" << endl;
	  continue;
      }
      psfBasis = gl.getBasis();
      psfBasis.setMu( psfBasis.getMu() + log(dx));
#ifdef DEBUGFDNTPSFEX
      cerr << "# PSF e and GL sigma: " << psfBasis.getS()
	   << " " << exp(psfBasis.getMu())
	   << endl;
#endif

      PSFInformation psfinfo(psfWCS, psfBasis);

      //
      // initialize galaxy ellipse value
      //
#ifdef DEBUGFDNTPSFEX
      cerr << "// (d) Starting ellipse of galaxy" << endl;
#endif

      cerr << "SExtractor ellipse r=" << r_pix*rescale << ", e="
	   << e1start << "," << e2start << endl;
      Ellipse sexE(e1start, e2start, log(r_pix*rescale), x_wc, y_wc);  // in wcs
      
      //
      // build FitExposure
      //
#ifdef DEBUGFDNTPSFEX
      cerr << "// (e) building FitExposure" << endl;
#endif
      FitExposure<>* fep;
#ifdef DEBUGFDNTPSFEX
      cerr << "// get the postage stamp" << endl;
#endif
      int xi0 = x_pix-stampSize/2;
      if (xi0<bounds.getXMin()) xi0=bounds.getXMin();
      if (xi0+stampSize+1>bounds.getXMax()) xi0 = bounds.getXMax()-stampSize-1;
      int yi0 = y_pix-stampSize/2;
      if (yi0<bounds.getYMin()) yi0=bounds.getYMin();
      if (yi0+stampSize+1>bounds.getYMax()) yi0 = bounds.getYMax()-stampSize-1;
#ifdef DEBUGFDNTPSFEX
      cerr << "// " << xi0 << " " << xi0+stampSize-1 << " " << yi0 << " " << yi0+stampSize-1 << endl;
      cerr << "// hopefully inside x=" << bounds.getXMin() << ".." << bounds.getXMax()
	   << ", y=" << bounds.getYMin() << ".." << bounds.getYMax() << endl;
#endif      
      Bounds<int> stamp(xi0,xi0+stampSize-1, yi0, yi0+stampSize-1);
      cerr << "set bounds" << endl;
      Image<float> scistamp = sci.subimage(stamp).duplicate();
      cerr << "sci bounds" << endl;
      Image<float> wtstamp  = wt.subimage(stamp).duplicate();
      cerr << "wt  bounds" << endl;
      Image<float> segstamp = seg.subimage(stamp).duplicate();
      cerr << "seg bounds" << endl;
      Image<double> xwstamp(stamp);
      Image<double> ywstamp(stamp);
      int segId=0;
      if (segIdCol>0)
	  segId = atoi(readvals[segIdCol-1].c_str());
      else {
	  bool badSegID=false;
	  for (int dx=0; dx<=1; dx++) {
	      for (int dy=0; dy<=1; dy++) {
		  cerr << int(x_pix)+dx << " " << int(y_pix)+dy << endl;
		  int dsegId=seg(int(x_pix)+dx,int(y_pix)+dy);
		  if (!segId && dsegId) {
		      segId=dsegId;
		  }
		  else if (dsegId && segId!=dsegId)
		  {
		      badSegID=true;
		  }
	      }
	  }

	  if (badSegID)
	  {
	      cerr << "# fail: could not get segmentation id for " << id << endl;
#ifdef CHECKPLOTS
	      ipsf.shift(1,1);
	      FITSImage<>::writeToFITS("check_"+id+"_psf_badseg.fits",ipsf);
	      wtstamp.shift(1,1);
	      FITSImage<>::writeToFITS("check_"+id+"_wt_badseg.fits",wtstamp);
	      scistamp.shift(1,1);
	      FITSImage<>::writeToFITS("check_"+id+"_sci_badseg.fits",scistamp);
	      segstamp.shift(1,1);
	      FITSImage<>::writeToFITS("check_"+id+"_seg_badseg.fits",segstamp);
#endif
	      nfail++;
	      nfail_seg++;
	      continue;
	  }
      }
      
      cerr << "got segid" << endl;
      {
	  try {
	      for (int iy=stamp.getYMin(); iy<=stamp.getYMax(); iy++)
	      {
		  for (int ix=stamp.getXMin(); ix<=stamp.getXMax(); ix++) {
		      double dxw,dyw;
		      fullmap->toWorld(ix,iy,dxw,dyw);
		      xwstamp(ix,iy)=dxw;
		      ywstamp(ix,iy)=dyw;
		      if (segstamp(ix,iy) && segstamp(ix,iy)!=segId) //it's another object!
		      {
			  cerr << "WARNING: objid " << id << " is fractured" << endl;
			  for (int iiy=max(stamp.getYMin(),iy-segmentationPadding);
			       iiy<=min(stamp.getYMax(),iy+segmentationPadding);
			       iiy++)
			      for (int iix=max(stamp.getXMin(),ix-segmentationPadding);
				   iix<=min(stamp.getXMax(),ix+segmentationPadding);
				   iix++)
				  wtstamp(iix,iiy)=0.;
		      }
		  }
	      }
	      // stack background with sky-subtracted single frames == photometric local background
	      // in stack flux scale
	      fep = new FitExposure<>(scistamp,
				      wtstamp,
				      xwstamp,
				      ywstamp, 0, sky+bg);
	}
	catch (...)
	{
	    cerr << "# fail: could not get postage stamp for " << id << endl;
	    nfail++;
	    nfail_post++;
	    delete fep; //delete map;
	    continue;
	}
      }

      //
      // check for bad pixels (there should be none in GREAT3!)
      //
#ifdef DEBUGFDNTPSFEX
      cerr << "// (f) checking bad pixels" << endl;
#endif
      double meanweight=0.;
      vector< Position<int> > bp = fep->badPixels(meanweight);
      cerr << "mean weight: " << meanweight << endl;
      if (bp.size() > maxBadPixels*stampSize*stampSize)
      {
	  cerr << "# fail: " << id << " has too many bad pixels" << endl;
#ifdef CHECKPLOTS
	  ipsf.shift(1,1);
	  FITSImage<>::writeToFITS("check_"+id+"_psf_badpix.fits",ipsf);
	  wtstamp.shift(1,1);
	  FITSImage<>::writeToFITS("check_"+id+"_wt_badpix.fits",wtstamp);
	  scistamp.shift(1,1);
	  FITSImage<>::writeToFITS("check_"+id+"_sci_badpix.fits",scistamp);
	  segstamp.shift(1,1);
	  FITSImage<>::writeToFITS("check_"+id+"_seg_badpix.fits",segstamp);
#endif
	  nfail++;
	  nfail_badpix++;
	  delete fep; // delete map;
	  continue;
      }
#ifdef DEBUGFDNTPSFEX
      // science image, from which we now subtract the GL model
      Image<float> glstamp = sci.subimage(stamp).duplicate();
#endif
      
      //
      // estimate "signal" of signal-to-noise (S/N) & interpolate over bad pixels
      //
      GLSimple<> gal(*fep, sexE, interpolationOrder);
      if (!gal.solve()) {
	  cerr << "# fail: GL fit failed for " << id << endl;
#ifdef CHECKPLOTS
	  ipsf.shift(1,1);
	  FITSImage<>::writeToFITS("check_"+id+"_psf_bpgl.fits",ipsf);
	  wtstamp.shift(1,1);
	  FITSImage<>::writeToFITS("check_"+id+"_wt_bpgl.fits",wtstamp);
	  scistamp.shift(1,1);
	  FITSImage<>::writeToFITS("check_"+id+"_sci_bpgl.fits",scistamp);
	  segstamp.shift(1,1);
	  FITSImage<>::writeToFITS("check_"+id+"_seg_bpgl.fits",segstamp);
#endif
	  nfail++;
	  nfail_bpgl++;
	  delete fep;
	  continue;
      }
      LVector bvec = gal.getB();
      double fluxModel = bvec.flux();   // this is the estimated signal
      Ellipse basis = gal.getBasis();
      

      // find half-light radius of observed galaxy
      double ee50obs = 1.0;
      // find the average rms... from the rms map!
      double imgRMS = 1.0;
      // calculate S/N based on half-light radius
      

      if (bp.size() > 0) // has bad pixels, but not too many to begin with: do GL interpolation
      {
#ifdef DEBUGFDNTPSFEX
	  cerr << "# trying to interpolate " << bp.size() << " bad pixels" << endl;
	  cerr << sexE << " " <<  interpolationOrder << endl;
#endif
	  double missingFlux=0.;
	  double scaleFactor = exp(basis.getMu());

	  for (vector< Position<int> >::iterator it=bp.begin(); it<bp.end(); it++) {
	      LVector psi(bvec.getOrder());
	      Position<double> xunit = basis.inv(Position<double>(xwstamp((*it).x,(*it).y),
								  ywstamp((*it).x,(*it).y)));
	      psi.fillBasis(xunit.x, xunit.y, scaleFactor);
	      scistamp((*it).x,(*it).y)=bvec.dot(psi);
	      wtstamp((*it).x,(*it).y)=meanweight;
	      // the interpolated pixel is assumed to have the mean weight of the good pixels;
	      // this makes sense because in the Fourier code homogeneous uncertainties are assumed
	      missingFlux += scistamp((*it).x,(*it).y);
	      scistamp((*it).x,(*it).y)+=sky+bg;
	  }
	  missingFlux *= rescale*rescale; // scale flux with coordinate system

#ifdef DEBUGFDNTPSFEX
	  for (int iy=stamp.getYMin(); iy<=stamp.getYMax(); iy++)
	  {
	      for (int ix=stamp.getXMin(); ix<=stamp.getXMax(); ix++) {
		  LVector psi(bvec.getOrder());
		  Position<double> xunit = basis.inv(Position<double>(xwstamp(ix,iy),
								      ywstamp(ix,iy)));
		  psi.fillBasis(xunit.x, xunit.y, scaleFactor);
		  glstamp(ix,iy) -= bvec.dot(psi);
	      }
	  }

	  cerr << "# fluxModel=" << fluxModel << endl;
	  cerr << "# scaleFactor=" << scaleFactor << endl;
	  cerr << "# missingFlux=" << missingFlux << endl;
#endif

	  if (missingFlux/fluxModel > maxBadFlux)
	  {
	      cerr << "# fail: bad pixels in " << id << " have too high flux fraction" << endl;
#ifdef CHECKPLOTS
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
#endif
	      nfail++;
	      nfail_bpff++;
	      delete fep;
	      continue;
	  }
      }

      //
      // run FDNT
      //
#ifdef DEBUGFDNTPSFEX      
      cerr << "// (g) running FDNT..." << flush;
#endif
      FDNT<> fd(*fep, psfinfo, sexE, order);
      fd.setMaskSigma(maskSigma);
      fd.GLAll();
      bool success = fd.prepare();
      if (success)
	  cerr << " success!" << endl;
      else
	  cerr << " failed." << endl;

      
#ifdef DEBUGFDNTPSFEX
      cerr << "// (h) evaluating FDNT" << endl;
#endif    
      Shear targetS = fd.getBasis().getS();
      cerr << targetS << endl;
      double prob;
      tmv::SymMatrix<double> covE(2,2);
      if (success) {
	  cerr << "making the actual measurement" << endl;
	  try {
	      targetS = fd.shape2(prob, covE);  // THIS IS THE ACTUAL MEASUREMENT!!
	  }
	  catch (...)
	  {
	      cerr << "an error occurred" << endl;
	  }
	  cerr << "made the actual measurement" << endl;
	  // ??? Make flag mask an input parameter ???
	  success = !(fd.getFlags() & (DidNotConverge + Singularity
				       + OutOfBounds + TooLarge + UnderSampled
				       + TooElliptical + GLFailure));
      }

      double eta1, eta2;
      targetS.getEta1Eta2(eta1,eta2);
      double egFix=0.;
      double sig1=0., sig2=0.;
      if (success) {
	  try {
	      egFix = fd.shrinkResponse(targetS);
	      cerr << "got shrink response" << endl;
	      sig1 = sqrt(covE(0,0));
	      sig2 = sqrt(covE(1,1));
	      se.add(targetS, egFix, sig1, sig2);
	      cerr << "good " << targetS << " " << egFix << " " << sig1 << " " << sig2 << endl;
	      ngood++;
	  }
	  catch (...)
	  {
	      cerr << "an error occured on previously successful measurement" << endl;
	      nfail++;
	      nfail_postmeas++;
	      success=false;
	  }
#ifdef CHECKPLOTS
	  ipsf.shift(1,1);
	  FITSImage<>::writeToFITS("check_"+id+"_psf_good.fits",ipsf);
	  wtstamp.shift(1,1);
	  FITSImage<>::writeToFITS("check_"+id+"_wt_good.fits",wtstamp);
	  scistamp.shift(1,1);
	  FITSImage<>::writeToFITS("check_"+id+"_sci_good.fits",scistamp);
	  segstamp.shift(1,1);
	  FITSImage<>::writeToFITS("check_"+id+"_seg_good.fits",segstamp);
	  Image<> filter1 = fd.drawFilter1(targetS, sexE);
	  Image<> filter2 = fd.drawFilter2(targetS, sexE);
	  filter1.shift(1,1);
	  FITSImage<>::writeToFITS("check_"+id+"_f1_good.fits",filter1);
	  filter2.shift(1,1);
	  FITSImage<>::writeToFITS("check_"+id+"_f2_good.fits",filter2);
	  if (bp.size()>0) {
	      glstamp.shift(1,1);
	      FITSImage<>::writeToFITS("check_"+id+"_glr_good.fits",glstamp);
	  }
#endif
      } else {
	  cerr << "# fail: some obscure thing in FDNT::prepare() happens for " << id
	       << ", flags " << fd.getFlags() << endl;
	stringstream flags; flags << fd.getFlags();	
#ifdef CHECKPLOTS
	ipsf.shift(1,1);
	FITSImage<>::writeToFITS("check_"+id+"_psf_fdnt_"+flags.str()+".fits",ipsf);
	wtstamp.shift(1,1);
	FITSImage<>::writeToFITS("check_"+id+"_wt_fdnt_"+flags.str()+".fits",wtstamp);
	scistamp.shift(1,1);
	FITSImage<>::writeToFITS("check_"+id+"_sci_fdnt_"+flags.str()+".fits",scistamp);
	segstamp.shift(1,1);
	FITSImage<>::writeToFITS("check_"+id+"_seg_fdnt_"+flags.str()+".fits",segstamp);
	if (bp.size()>0) {
	glstamp.shift(1,1);
	FITSImage<>::writeToFITS("check_"+id+"_glr_fdnt_"+flags.str()+".fits",glstamp);
	}
#endif
	
	nfail++;
	nfail_fdnt++;
      }
      double mu = fd.getBasis().getMu();
      cout << id 
	   << " " << x_pix
	   << " " << y_pix 
	   << " " << obs_SN
	   << " " << eta1 
	   << " " << eta2
	   << " " << sig1
	   << " " << sig2
	   << " " << mu
	   << " " << egFix
	   << " " << fd.getFlags()
	   << " " << success
	   << " " << fwd
	   << endl;

      delete fep; //delete map;
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
    cout << "# Reasons for failure:\n# " << nfail_post << "\t couldn't get postage stamp\n# " << nfail_badpix << "\t had too many bad pixels\n# " 
         << nfail_bpgl << "\t couldn't GL-interpolate bad pixels\n# " << nfail_bpff << "\t had too much flux in bad pixels\n# " 
         << nfail_postmeas << "\t failed after measurement\n# "  
         << nfail_fdnt << "\t had some error in FDNT\n# " << nfail_seg << "\t had ambiguous segmentation information\n"<< endl;
  } catch (std::runtime_error &m) {
    quit(m,1);
  }
}
