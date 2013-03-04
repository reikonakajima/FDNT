// $Id: FDNTGreat.cpp,v 1.6 2011/07/15 17:46:17 garyb Exp $

// Run FDNT on GREAT images with SExtractor input files

#include "FDNT.h"
#include "StringStuff.h"
#include "Pset.h"
#include "FITSImage.h"
#include <fstream>
#include "ShearEstimator.h"

using namespace laguerre;
using namespace sbp;

int
main(int argc, 
     char *argv[])
{
  // Parameters to set:
  // columns for necessary input:
  int idCol;
  int xCol;
  int yCol;
  int magCol;
  int rCol;
  int aCol;
  int bCol;
  int paCol;

  // The fit:
  int order;
  double maskSigma;

  // Input files
  string fitsName;
  string catName;
  string psfName;

  // Data properties
  double sky;
  double rn;
  double gain;
  int stampSize;

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
			   "Image FITS file", "set0001.fit");
      parameters.addMember("catName",&catName, def,
			   "Input SExtractor catalog", "set0001.scat");
      parameters.addMember("psfName",&psfName, def,
			   "Input GL PSF file", "set0001.bvec");
      parameters.addMember("stampSize",&stampSize, def,
			   "Pixel size of img postage stamps (0=no bound)", 64);
      parameters.addMemberNoValue("DATA PROPERTIES:",0,
				  "Input data");
      parameters.addMember("sky",&sky, def,
			   "sky level", 0.);
      parameters.addMember("rn",&rn, def | low,
			   "Read noise, i.e. RMS at 0 ADU", 1000., 0.);
      parameters.addMember("gain",&gain, def | low,
			   "Gain, e/ADU for Poissson, 0 for no Poisson", 1., 0.);
      parameters.addMemberNoValue("FDNT:",0,
				  "Fitting parameters:");
      parameters.addMember("order",&order, def,
			   "FDNT GL-fit order (0 for FFT)", 0);
      parameters.addMember("maskSigma",&maskSigma, def,
			   "GL and FDNT mask sigma (-1 auto)", 4.);
      parameters.addMemberNoValue("CATALOG TRANSLATION:",0,
				  "Columns holding input data");
      parameters.addMember("idCol",&idCol, def | low,
			   "Object ID", 1, 1);
      parameters.addMember("xCol",&xCol, def | low,
			   "x centroid", 2, 1);
      parameters.addMember("yCol",&yCol, def | low,
			   "y centroid", 3, 1);
      parameters.addMember("magCol",&magCol, def | low,
			   "Magnitude", 4, 1);
      parameters.addMember("rCol",&rCol, def | low,
			   "Half-light radius", 8, 1);
      parameters.addMember("aCol",&aCol, def | low,
			   "Major axis", 9, 1);
      parameters.addMember("bCol",&bCol, def | low,
			   "Minor axis", 10, 1);
      parameters.addMember("paCol",&paCol, def | low,
			   "Position angle", 11, 1);
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

    cout << "# " << stringstuff::taggedCommandLine(argc, argv) << endl;
    parameters.dump(cout);

    // Read the PSF
    SBProfile* psf=0;
    double e1psf, e2psf, mupsf;
    if (stringstuff::nocaseEqual(psfName,"moffat")) {
      //**** Make Moffat preloaded:
      double g1 = -0.03;  // negative is DEBUG, sign flips in some conventions!  (e.g., Great08)
      double g2 = -0.07;
      double g = hypot(g1,g2);
      double psfe=tanh(2*atanh(g));
      SBMoffat sbm(3.0, 2.);
      sbm.setFWHM(2.85/sqrt(1-g*g));
      SBProfile* sb2 = sbm.shear(g1*psfe/g, g2*psfe/g);
      SBBox box(1.);
      psf = new SBConvolve(*sb2, box);
      e1psf = -0.056;    // negative is DEBUG, sign flips in some conventions!  (e.g., Great08)
      e2psf = -0.12;
      mupsf = 0.344;    // probably need to calculate this properly... or not?
    } else {
      LVector bPSF;
      ifstream cpsf(psfName.c_str());
      if (!cpsf) {
	cerr << "Could not open PSF file " << psfName << endl;
	exit(1);
      }
      string buffer;
      stringstuff::getlineNoComment(cpsf, buffer);
      istringstream iss(buffer);
      if (!(iss >> e1psf >> e2psf >> mupsf)) {
	cerr << "Error on PSF ellipse: " << buffer << endl;
	exit(1);
      }
      cpsf >> bPSF;

      SBLaguerre sbl(bPSF, exp(mupsf));
      psf = sbl.shear(e1psf, e2psf);
    }

    PSFInformation psfinfo(*psf, Ellipse(e1psf, e2psf, mupsf));

    // Open image
    Image<> sci;
    FITSImage<> fits(fitsName);
    sci = fits.extract();
      
    // create coordinate maps and weight map
    Bounds<int> b=sci.getBounds();
    Image<> wt(b, pow(rn,-2.));
    Image<double> xw(b);
    Image<double> yw(b);
    for (int iy=b.getYMin(); iy<=b.getYMax(); iy++)
      for (int ix=b.getXMin(); ix<=b.getXMax(); ix++) {
	xw(ix,iy) = ix;
	yw(ix,iy) = iy;
	// ??? Following not quite right for gain !=1 :
	if (gain>0.) wt(ix,iy) = 1./ (rn*rn + MAX(sci(ix,iy), (float) 0.) / gain);
      }

    // initialize sums
    UnweightedShearEstimator se;
    int ngood=0;
    int nfail=0;

    // Loop through objects
    int nread=MAX(idCol, xCol);
    nread = MAX(nread, yCol);
    nread = MAX(nread, magCol);
    nread = MAX(nread, rCol);
    nread = MAX(nread, aCol);
    nread = MAX(nread, bCol);
    nread = MAX(nread, paCol);
    vector<string> readvals(nread);
    ifstream ccat(catName.c_str());
    if (!ccat) {
      cerr << "Error opening catalog file " << catName << endl;
      exit(1);
    }
    string buffer;
    tmv::SymMatrix<double> covE(2);

    while (stringstuff::getlineNoComment(ccat, buffer)) {
      // Acquire info
      istringstream iss(buffer);
      for (int i=0; i<nread; i++) iss >> readvals[i];
      if (!iss) {
	cerr << "Bad catalog input line: " << buffer;
	exit(1);
      }
      string id = readvals[idCol-1];
      double x0 = atof(readvals[xCol-1].c_str());
      double y0 = atof(readvals[yCol-1].c_str());
      double mag = atof(readvals[magCol-1].c_str());
      double r = atof(readvals[rCol-1].c_str());
      double a = atof(readvals[aCol-1].c_str());
      double b = atof(readvals[bCol-1].c_str());
      double pa = atof(readvals[paCol-1].c_str());

      // Starting ellipse
      double e1start = (a*a-b*b)/(a*a+b*b);
      double e2start = e1start * sin(2*pa*PI/180.);
      e1start *= cos(2*pa*PI/180.);

      Ellipse sexE(e1start, e2start, log(r), x0, y0);

      FitExposure<>* fep;
      if (stampSize>0) {
        /* old code, does not make sense - OM
	int x1 = stampSize*static_cast<int> (floor(x0/stampSize))+1;
	int y1 = stampSize*static_cast<int> (floor(y0/stampSize))+1;
	*/
	Bounds<int> sci_bd=sci.getBounds();
	int x1 = x0 - static_cast<int> (floor(stampSize/2));
	if ( x1 < 1 ) {
	  x1 = 1;
	}
	if ( x1+stampSize-1 > sci_bd.getXMax() ) {
	  x1 = sci_bd.getXMax() - (stampSize-1);
	}
	// will fail if stampSize > ( b.getXMax - b.getXMin )  !!!
	int y1 = y0 - static_cast<int> (floor(stampSize/2));
	if ( y1 < 1 ) {
	  y1 = 1;
	}
	if ( y1+stampSize-1 > sci_bd.getYMax() ) {
	  y1 = sci_bd.getYMax() - (stampSize-1);
	}
	// will fail if stampSize > ( b.getXMax - b.getXMin )  !!!
	Bounds<int> stamp(x1,x1+stampSize-1, y1, y1+stampSize-1);
	fep = new FitExposure<>(sci.subimage(stamp), 
				wt.subimage(stamp), 
				xw.subimage(stamp), 
				yw.subimage(stamp), 0, sky);
      } else {
	fep = new FitExposure<>(sci, wt, xw, yw, 0, sky);
      }
      FDNT<> fd(*fep, psfinfo, sexE, order);
      fd.setMaskSigma(maskSigma);
      fd.GLAll();
      bool success = fd.prepare();

      Shear targetS = fd.getBasis().getS();
      double prob;
      tmv::SymMatrix<double> covE(2,2);
      if (success) {
	try {
	  targetS = fd.shape2(prob, covE);  // THIS IS THE ACTUAL MEASUREMENT!!
	  // ??? Make flag mask an input parameter ???
	  success = !(fd.getFlags() & (DidNotConverge + Singularity
				       + OutOfBounds + TooLarge + UnderSampled
				       + TooElliptical + GLFailure));
	} catch (std::runtime_error &m) {
	  success = false;
	  fd.setFlag(Singularity);
	}
      /**** old stuff:
	fd.doglegNull(0.01);
	targetS = fd.getBasis().getS();
	fd.resizeGalaxyWt(targetS);
	success = fd.doglegNull(0.01);
	targetS = fd.getBasis().getS();
	fd.resizeGalaxyWt(targetS);

	double prob;
	targetS = fd.shape(prob, covE);
	success = !(fd.getFlags() & (DidNotConverge + Singularity
				     + OutOfBounds + TooLarge
				     + TooElliptical + GLFailure));
      ***/
      }

      double eta1, eta2;
      targetS.getEta1Eta2(eta1,eta2);
      double egFix=0.;
      double sig1=0., sig2=0.;
      if (success) {
	egFix = fd.shrinkResponse(targetS);
	sig1 = sqrt(covE(0,0));
	sig2 = sqrt(covE(1,1));
	se.add(targetS, egFix, sig1, sig2);
	ngood++;
      } else {
	nfail++;
      }
      double mu = fd.getBasis().getMu();
      cout << id 
	   << " " << x0 
	   << " " << y0 
	   << " " << mag 
	   << " " << eta1 
	   << " " << eta2
	   << " " << sig1
	   << " " << sig2
	   << " " << mu
	   << " " << egFix
	   << " " << fd.getFlags()
	   << endl;

      delete fep;
    }

    // Print out mean shear estimate
    Shear S(se);
    double g1, g2, sig1, sig2;
    S.getG1G2(g1,g2);
    se.sigmaE(sig1,sig2, false);

    cout << fixed << setprecision(6);
    // Approximate the reduced-shear error as 1/2 of the distortion error
    cout << "#Means: " << g1 << " +- " << sig1/2
	 << " " << g2 << " +- " << sig2/2
	 << endl;
    cout << "#Failed / good " << nfail << " / " << ngood << endl;
  } catch (std::runtime_error &m) {
    quit(m,1);
  }
}
