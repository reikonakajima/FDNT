// $Id: FDNTEllipse.cpp,v 1.6 2011/07/18 15:36:07 dgru Exp $

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
  // The ring:
  int nDither;
  int realizationsPerDither;

  // The fit:
  int order;
  double fixedSigma;
  bool fixSigma;
  int iFixCentroid;
  bool fixCentroid;

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

      parameters.addMemberNoValue("TEST:",0,
				  "Test specifications");
      parameters.addMember("nDither",&nDither, def | low,
			   "Number of dithers", 10, 1);
      parameters.addMember("nNoise",&realizationsPerDither, def | low,
			   "Number of noise realizations per dither", 100, 1);
      parameters.addMemberNoValue("SCENE:",0,
				  "Galaxy and PSF to test");
      parameters.addMember("sbGalaxy",&sbGalaxy, def,
			   "SBParse string for galaxy", "exp");
      parameters.addMember("sig",&significance, def | low,
			   "Obs. detection S/N within half-light circle", 100., 1.);
      parameters.addMember("rhalf",&rhalf, def,
			   "Post-shear gal. half-light rad. (<=0 to use input size)", -1.);
      parameters.addMember("e1",&e1in, def | lowopen | upopen, 
			   "Applied ellipticity + component", 0.3, -1., 1.);
      parameters.addMember("e2",&e2in, def | lowopen | upopen, 
			   "Applied ellipticity x component", 0.0, -1., 1.);
      parameters.addMember("sbPSF",&sbPSF, def,
			   "SBParse string for PSF", "gauss * box");
      parameters.addMember("rhalfPSF",&rhalfPSF, def,
			   "PSF half-light radius (<=0 to use input size)", -1.);
      parameters.addMemberNoValue("FDNT:",0,
				  "Fitting parameters:");
      parameters.addMember("order",&order, def,
			   "FDNT GL-fit order (0 for FFT)", 0);
      parameters.addMember("fixedSigma",&fixedSigma, def,
			   "Fixed galaxy weight function sigma (<=0 to fit for value)", 0.);
      parameters.addMember("fixCentroid",&iFixCentroid, def,
			   "Fix centroid at correct value?", 0);
      parameters.addMember("catfile", &catfile, def,
			   "Output catalog file", "junk.cat");
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

    // Will we be letting galaxy sigma vary?
    fixSigma = (fixedSigma > 0.);
    fixCentroid = (iFixCentroid <= 0.);

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

      /*
      // save the object as a FITS file
      clean.shift(1,1);
      FITSImage<>::writeToFITS("clean.fits", clean);
      */
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

    ofs << "# eta1  eta2  +- err1  err2  r12  Sigma  EG   Log(L)" << endl;

    for (int iDither=0; iDither<nDither; iDither++) {
      double xshift = u-0.5;
      double yshift = u-0.5;
      SBProfile* gshift = galaxy->shift(xshift,yshift);
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

      // Noise realization loop:
      int iNoise=0;
      int failNoise=0;	// number of failures at this step
      while (iNoise<realizationsPerDither) {
	Image<> sci = noiseless.duplicate();
	for (int iy=sci.YMin(); iy<=sci.YMax(); iy++)
	  for (int ix=sci.XMin(); ix<=sci.XMax(); ix++) 
	    sci(ix,iy) += g*imgRMS; //******???*/ *0.01;
	
	FitExposure<> fe(sci, wt, xw, yw, 0);
	if (fixCentroid) {
	  // Give exact center, assuming PSF and galaxy both
	  // had centers at origin before the shift:
	  startBasis.setX0(Position<double>(xshift, yshift));
	}
	FDNT<> fd(fe, psfinfo, startBasis, order);
	if (fixCentroid) fd.setCentering(false);
	fd.GLAll();  // Measure observed size/centroid with GL fitting
	fd.setMaskSigma(6.);     // Decide how much of image to be used for FDNT
	bool success = fd.prepare();  // Extract Fourier info for FDNT
	//	if (success) success = fd.doglegNull(0.01);
	if (fixSigma) {
	  fd.setIntrinsicSize(fixedSigma);
	  fd.setSizing(false);
	}
	Shear targetS = fd.getBasis().getS();  
	double prob;
	tmv::SymMatrix<double> covE(2,2);
	if (success) {
	  targetS = fd.shape2(prob, covE);  // THIS IS THE ACTUAL MEASUREMENT!!
	  // ??? Make flag mask an input parameter ???
	  success = !(fd.getFlags() & (DidNotConverge + Singularity
				       + OutOfBounds + TooLarge + UnderSampled
				       + TooElliptical + GLFailure));
	}
	if (!success) {
	  failNoise++;
	  nFail++;
	  if (failNoise > 10 && failNoise > 4*realizationsPerDither) {
	    // ?? Check limits above - quit if ~>50% of trials failing 
	    cerr << "Too many failures at one dither" << endl;
	    exit(1);
	  }
	  continue;
	}

	// Success: accumulate statistics and output data
	double sigma = exp(fd.getBasis().getMu());
	double egFix = fd.shrinkResponse(targetS);
	double eta1, eta2;
	targetS.getEta1Eta2(eta1,eta2);

	ofs << setw(8) << eta1
	    << " " << setw(8) << eta2
	    << " " << setw(8) << sqrt(covE(0,0))
	    << " " << setw(8) << sqrt(covE(1,1))
	    << " " << setw(8) << covE(1,0)/sqrt(covE(0,0)*covE(1,1))
	    << " " << setw(8) << sigma
	    << " " << setw(8) << egFix
	    << " " << setw(5) << prob
	    << " " << setw(4) << fd.getFlags()
	    << endl;

	/**/ofs << "# " << fd.getEvaluationCount()
		<< " evals in " << fd.getETrialCount()
		<< " trials, density " << fd.getSampleDensity()
		<< " totprob " << fd.getTotalProbability()
		<< endl;

	iNoise++;
	nGood++;
	eta1 -= eta1in;
	eta2 -= eta2in;
	sum1 += eta1; sum2+=eta2;
	sum11 += eta1*eta1; sum12+=eta1*eta2; sum22+=eta2*eta2;
	sum111 += pow(eta1,3.); sum222+=pow(eta2,3.);
	sum1111 += pow(eta1,4.); sum2222+=pow(eta2,4.);
	
	sumEG += egFix; sumEG2+=egFix*egFix;
	sumLL += prob; sumLL2 += prob*prob;
	sumSigma += sigma; sumSigma2+=sigma*sigma;
	sumVar1 += covE(0,0); sumVar2+= covE(1,1); sumCov12+= covE(0,1);

      } // end realization loop
    } // end dither loop

    // Print out statistics 
    cout << "N " << nGood << endl;
    cout  << "Failure rate: " << (1.*nFail) / (nFail+nGood) << endl;
    double bias1 = sum1 / nGood;
    double bias2 = sum2 / nGood;
    double var1 = sum11/nGood - bias1*bias1;
    double var2 = sum22/nGood - bias2*bias2;
    double cov12 = sum12/nGood - bias1*bias2;
    cout << "Bias " << setw(8) << bias1
	 << " "  << setw(8) << bias2
	 << " +- " << setw(8) << sqrt(var1/nGood)
	 << " " << setw(8) << sqrt(var2/nGood)
	 << endl;

    cout << "RMS, r: " << setw(8) << sqrt(var1)
	 << " " << setw(8) << sqrt(var2)
	 << " " << setw(8) << cov12/sqrt(var1*var2)
	 << endl;
    cout << "Cov elements " << setw(8) << var1
	 << " " << setw(8) << var2
	 << " " << setw(8) << cov12
	 << endl;
    cout << "Estimated " << setw(8) << sumVar1/nGood
	 << " " << setw(8) << sumVar2/nGood
	 << " " << setw(8) << sumCov12/nGood
	 << endl;

    cout << "Skewness " << setw(5) << (sum111/nGood - 3*bias1*var1 - pow(bias1,3.))*pow(var1, -1.5) 
	 << " " << setw(5) << (sum222/nGood - 3*bias2*var2 - pow(bias2,3.))*pow(var2, -1.5) 
	 << endl;

    cout << "Kurtosis " << setw(5) 
	 << (sum1111/nGood - 4*sum111*bias1/nGood + 6*var1*bias1*bias1 + 3*pow(bias1, 4.)) * pow(var1,-2.) - 3.
	 << " " << setw(5) 
	 << (sum2222/nGood - 4*sum222*bias2/nGood + 6*var2*bias2*bias2 + 3*pow(bias2, 4.)) * pow(var2,-2.) - 3.
	 << endl;

    cout << "Log(L) " << sumLL/nGood 
	 << " +- " << sqrt( sumLL2/nGood - pow(sumLL/nGood, 2.))
	 << endl;

    cout << "EG " << sumEG/nGood 
	 << " +- " << sqrt( sumEG2/nGood - pow(sumEG/nGood, 2.))
	 << endl;

    cout << "sigma " << sumSigma/nGood
	 << " +- " << sqrt( sumSigma2/nGood - pow(sumSigma/nGood, 2.))
	 << endl;

    delete galaxy;
    delete psfraw;
  } catch (std::runtime_error &m) {
    quit(m,1);
  }
}

