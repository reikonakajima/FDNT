// $Id: FDNTEllipse.cpp,v 1.1 2010/04/28 18:31:32 garyb Exp $

// Run many trials of FDNT measure of a single
// galaxy shape

#include "FDNT.h"
#include "StringStuff.h"
#include "Pset.h"
#include "FITSImage.h"
#include "Random.h"
#include "SBParse.h"
#include "EnclosedFluxRadius.h"

using namespace laguerre;
using namespace sbp;
using namespace ran;

int
main(int argc, 
     char *argv[])
{
  // Parameters to set:
  // The grid:
  double estep;
  double e1min, e1max, e2min, e2max;

  // The fit:
  int order;

  // The galaxy - specified by an SBParse string
  string sbGalaxy;
  double significance;
  // Scale to specified half-light radius if desired:
  double rhalf;

  // The PSF - from an SBParse string
  string sbPSF;
  // Scale to specified half-light radius if desired:
  double rhalfPSF;

  // Applied shear to basic galaxy
  double e1in;
  double e2in;

  string catfile;

  Pset parameters;
  {
      const int def=PsetMember::hasDefault;
      const int low=PsetMember::hasLowerBound;
      const int up=PsetMember::hasUpperBound;
      const int lowopen = low | PsetMember::openLowerBound;
      const int upopen = up | PsetMember::openUpperBound;

      parameters.addMemberNoValue("GRID:",0,
				  "grid of ellipticities to probe:");
      parameters.addMember("e1min",&e1min, def | low | up,
			   "Minimum e1", -1., -1., 1.);
      parameters.addMember("e1max",&e1max, def | low | up,
			   "Maximum e1", -1., -1., 1.);
      parameters.addMember("e2min",&e2min, def | low | up,
			   "Minimum e2", -1., -1., 1.);
      parameters.addMember("e2max",&e2max, def | low | up,
			   "Maximum e2", -1., -1., 1.);
      parameters.addMember("estep",&estep, def | lowopen,
			   "e grid step", 0.1, 0.);
      parameters.addMember("sbGalaxy",&sbGalaxy, def,
			   "SBParse string for galaxy", "exp");
      parameters.addMember("sig",&significance, def | low,
			   "Obs. detection S/N within half-light circle", 1e4, 1.);
      parameters.addMember("rhalf",&rhalf, def,
			   "Post-shear gal. half-light rad. (<=0 to use input size)", -1.);
      parameters.addMember("e1",&e1in, def | lowopen | upopen, 
			   "Applied ellipticity + component", 0.3, -1., 1.);
      parameters.addMember("e2",&e2in, def | lowopen | upopen, 
			   "Applied ellipticity x component", 0.0, -1., 1.);
      parameters.addMember("sbPSF",&sbPSF, def,
			   "SBParse string for PSF", "gauss");
      parameters.addMember("rhalfPSF",&rhalfPSF, def,
			   "PSF half-light radius (<=0 to use input size)", -1.);
      parameters.addMemberNoValue("FDNT:",0,
				  "Fitting parameters:");
      parameters.addMember("order",&order, def,
			   "FDNT GL-fit order (0 for FFT)", 0);
      parameters.addMember("catfile", &catfile, def,
			   "Output catalog file", "/dev/null");
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

    // Open output file for catalog
    ofstream ofs(catfile.c_str());
    if (!ofs) {
	cerr << "Can't open output catalog file " << catfile << endl;
	exit(1);
    }

    cout << "# " << stringstuff::taggedCommandLine(argc, argv) << endl;
    parameters.dump(cout);
    ofs << "# " << stringstuff::taggedCommandLine(argc, argv) << endl;
    parameters.dump(ofs);


    int nx = static_cast<int> ( ceil( (e1max-e1min)/estep+0.5) );
    int ny = static_cast<int> ( ceil( (e2max-e2min)/estep+0.5) );
    Bounds<int> b(1,nx,1,ny);
    Image<> probImg(b, 0.);
    Image<> sigImg(b, 0.);
    Image<> detImg(b, 0.);

    // Start up our random number generators
    UniformDeviate u;
    GaussianDeviate g;


    // First create the galaxy
    SBProfile* unlensed=SBParse(sbGalaxy);
    // ..and make it elliptical
    SBProfile* galaxy = unlensed->shear(e1in,e2in);
    delete unlensed; unlensed=0;

    // Measure half-light radius
    double ee50g = EnclosedFluxRadius(*galaxy);
    // Dilate if desired:
    if (rhalf > 0.) {
      cout << "# Galaxy half-light radius " << ee50g << endl;
      double factor = rhalf / ee50g;
      Ellipse e(0.,0.,log(factor));
      SBProfile* sb2 = galaxy->distort(e);
      delete galaxy;
      galaxy = sb2;
      ee50g = rhalf;
    }
    cout << "# Galaxy EE50: " << ee50g << endl;

    // Create PSF
    SBProfile* psfraw=SBParse(sbPSF);
    double ee50psf = EnclosedFluxRadius(*psfraw);
    // Dilate if desired:
    if (rhalfPSF > 0.) {
      cout << "# Raw PSF half-light radius " << ee50psf << endl;
      double factor = rhalfPSF / ee50psf;
      Ellipse e(0.,0.,log(factor));
      SBProfile* sb2 = psfraw->distort(e);
      delete psfraw;
      psfraw = sb2;
      ee50psf = rhalfPSF;
    }
    psfraw->setFlux(1.);	// Always unit normalized
    cout << "# PSF EE50: " << ee50psf << endl;

    // Draw a version of the object at zero rotation to decide noise
    // level that makes desired significance
    double gflux=-1.;
    double imgRMS = -1.;
    double estSize = hypot( ee50g, ee50psf);
    Ellipse cleanBasis(0., 0., log(estSize));
    {
      SBConvolve obs(*galaxy, *psfraw);
      double ee50obs = EnclosedFluxRadius(obs);
      cout << "# Observed EE50: " << ee50obs << endl;
      gflux = obs.getFlux();
      cout << "# Galaxy flux: " << gflux << endl;
      // Set noise level from circular-aperture formula:
      imgRMS = 0.5*gflux / (sqrt(PI)*ee50obs*significance);
      cout << "# Pixel RMS: " << imgRMS << endl;
      
      // Measure GL size & significance
      // Do we need to oversample it?
      double dx = (ee50obs > 3.) ? 1. : 0.5;
      Image<> clean = obs.draw(dx);
      clean.getHdrValue("DX",dx);
      cleanBasis.setMu( cleanBasis.getMu() - log(dx));

      Image<> wt(clean.getBounds());
      // Set wt for very high S/N:
      const double trialSN=1e5;
      wt = pow(imgRMS*significance/(dx*trialSN), -2.);
      GLSimple<> gl(clean, wt, cleanBasis, 4);
      // Not in general:
      //  gl.setCentering(false);
      if (!gl.solve()) {
	cerr << "Failed measuring clean galaxy, flags " << gl.getFlags() << endl;
	clean.shift(1,1);
	FITSImage<>::writeToFITS("clean.fits", clean);
	exit(1);
      }
      cleanBasis = gl.getBasis();
      cleanBasis.setMu( cleanBasis.getMu() + log(dx));
      cout << "# Observed GL sigma: " << exp(cleanBasis.getMu()) << endl;
      double f, varf;
      gl.b00(f, varf);
      cout << "# Observed GL S/N: " << (f / sqrt(varf)) * (significance / trialSN) << endl;

      // Calculate optimal detection significance
      double sigsum=0.;
      double fluxsum=0.;
      for (int iy=clean.YMin(); iy<=clean.YMax(); iy++)
	for (int ix=clean.XMin(); ix<=clean.XMax(); ix++) {
	  sigsum += clean(ix,iy)*clean(ix,iy);
	  fluxsum += clean(ix,iy);
	}
      cout << "# Observed ideal S/N: " << sqrt(sigsum)*dx/imgRMS << endl;

      // save the object as a FITS file
      clean.shift(1,1);
      FITSImage<>::writeToFITS("clean.fits", clean);
    }

    // Make a PSFInfo structure - do a GL fit to get its basis
    Ellipse psfBasis(0., 0., log(ee50psf/1.17));
    {
      double dx = MAX(ee50psf, 1.) / 4.;
      dx = MIN(1., dx);
      Image<> ipsf = psfraw->draw(dx);

      // Measure GL size & significance
      Image<> wt(ipsf.getBounds());
      wt = pow(1e-5/dx, -2.);
      psfBasis.setMu( psfBasis.getMu() - log(dx));
      GLSimple<> gl(ipsf, wt, psfBasis, 4);
      if (!gl.solve()) {
	cerr << "Failed measuring psf, flags " << gl.getFlags() << endl;
	ipsf.shift(1,1);
	FITSImage<>::writeToFITS("ipsf.fits", ipsf);
	exit(1);
      }
      psfBasis = gl.getBasis();
      psfBasis.setMu( psfBasis.getMu() + log(dx));
      cout << "# PSF e and GL sigma: " << psfBasis.getS()
	   << " " << exp(psfBasis.getMu())
	   << endl;
    }
    PSFInformation psfinfo(*psfraw, psfBasis);

    Ellipse startBasis = cleanBasis;

    // Collect stats:
    long int nFail = 0;	// Number of failures
    long int nGood=0;
    double sum1=0., sum2=0., sum11=0., sum12=0.;
    double sum22=0., sum111=0., sum222=0., sum1111=0., sum2222=0.;
    double sumEG = 0., sumEG2=0.;
    double sumLL = 0., sumLL2 = 0.;
    double sumVar1=0., sumVar2=0., sumCov12=0.;
    double sumSigma=0., sumSigma2=0.;
    double eta1in, eta2in;
    {
      Shear shearin(e1in,e2in);
      shearin.getEta1Eta2(eta1in, eta2in);
    }

    SBProfile* gshift = galaxy->shift(u-0.5, u-0.5);
    SBConvolve final(*psfraw, *gshift);

    // Draw at 1 unit per pixel
    Image<> noiseless = final.draw(1.);
    // Don't need gshift any more:
    delete gshift; gshift=0;

    // Make our weight and coordinate images:
    Image<> wt(noiseless.getBounds(), pow(imgRMS, -2.));
    Image<double> xw(noiseless.getBounds());
    Image<double> yw(noiseless.getBounds());
    for (int iy=xw.YMin(); iy<=xw.YMax(); iy++)
      for (int ix=xw.XMin(); ix<=xw.XMax(); ix++) {
	xw(ix,iy) = ix;
	yw(ix,iy) = iy;
      }

    Image<> sci = noiseless.duplicate();
    for (int iy=sci.YMin(); iy<=sci.YMax(); iy++)
      for (int ix=sci.XMin(); ix<=sci.XMax(); ix++) 
	sci(ix,iy) += g*imgRMS;
	
    FitExposure<> fe(sci, wt, xw, yw, 0);
    FDNT<> fd(fe, psfinfo, startBasis, order);
    fd.GLAll();
    fd.setMaskSigma(6.);
    fd.prepare();

    for (int ix=1; ix<=nx; ix++) {
      for (int iy=1; iy<=ny; iy++) {
	double e1=e1min + estep*(ix-1);
	double e2=e2min + estep*(iy-1);
	if (hypot(e1,e2)>=1.) continue;
	Shear targetS(e1,e2);
		      
	tmv::SymMatrix<double> covE(2,2);
	double prob = fd.logProbability(targetS, covE);
	// Success: accumulate statistics and output data
	double sigma = exp(fd.getBasis().getMu());
	double det = covE.Det();
	probImg(ix,iy) = prob;
	sigImg(ix,iy) = sigma;
	detImg(ix,iy) = det;

	ofs << setw(8) << e1
	    << " " << setw(8) << e2
	    << " " << setw(8) << prob
	    << " " << setw(8) << sigma
	    << " " << setw(8) << det
	    << endl;

      } // end iy loop
    } // end ix loop
    FITSImage<>::writeToFITS("prob.fits",probImg);
    FITSImage<>::writeToFITS("sig.fits",sigImg);
    FITSImage<>::writeToFITS("det.fits",detImg);

    delete galaxy;
    delete psfraw;
  } catch (std::runtime_error &m) {
    quit(m,1);
  }
}

