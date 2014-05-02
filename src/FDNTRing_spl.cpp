// $Id: FDNTRing.cpp,v 1.10 2011/10/20 02:37:08 garyb Exp $

// Ring test with FDNT code

#include "FDNT.h"
#include "StringStuff.h"
#include "Pset.h"
#include "FITSImage.h"
#include "Random.h"
#include "ShearEstimator.h"
#include "SBParse.h"
#include "EnclosedFluxRadius.h"
#include "AstronomicalConstants.h"

using namespace laguerre;
using namespace sbp;
using namespace ran;

int
main(int argc, 
     char *argv[])
{
  // Parameters to set:
  // The ring:
  int nTheta;  // number of evenly spaced directions
  int iTheta;  // the ith direction
  int dithersPerTheta;
  int realizationsPerDither;

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

  // Applied shear
  double dist1;
  double dist2;

  // more measurement outputs?
  // 0: no extra outputs
  // 1: individual shapes (relevant for noisy measurements)
  int    outputFlags;
  string outputFilename;

  Pset parameters;
  {
      const int def=PsetMember::hasDefault;
      const int low=PsetMember::hasLowerBound;
      const int up=PsetMember::hasUpperBound;
      const int lowopen = low | PsetMember::openLowerBound;
      const int upopen = up | PsetMember::openUpperBound;

      parameters.addMemberNoValue("RING:",0,
				  "Ring test specifications");
      parameters.addMember("nTheta",&nTheta, def | low,
			   "Number of orientations", 32, 4);
      parameters.addMember("iTheta",&iTheta, def | low,
			   "The ith orientation", 0, 0);
      parameters.addMember("nDither",&dithersPerTheta, def | low,
			   "Number of dithers per orientation", 3, 1);
      parameters.addMember("nNoise",&realizationsPerDither, def | low,
			   "Number of noise realizations per dither", 10, 1);
      parameters.addMember("dist1",&dist1, def | lowopen | upopen, 
			   "Applied distortion + component", 0.01, -1., 1.);
      parameters.addMember("dist2",&dist2, def | lowopen | upopen, 
			   "Applied distortion x component", 0.0, -1., 1.);
      parameters.addMember("sbGalaxy",&sbGalaxy, def,
			   "SBParse string for galaxy", "exp");
      parameters.addMember("sig",&significance, def | low,
			   "Obs. detection S/N within half-light circle", 1e4, 1.);
      parameters.addMember("rhalf",&rhalf, def,
			   "Galaxy half-light radius (<=0 to use input size)", -1.);
      parameters.addMember("sbPSF",&sbPSF, def,
			   "SBParse string for PSF", "gauss");
      parameters.addMember("rhalfPSF",&rhalfPSF, def,
			   "PSF half-light radius (<=0 to use input size)", -1.);
      parameters.addMemberNoValue("FDNT:",0,
				  "Fitting parameters:");
      parameters.addMember("order",&order, def,
			   "FDNT GL-fit order (0 for FFT)", 0);
      parameters.addMember("outputFlags", &outputFlags, def | low,
			   "Do we want outputs for individual measurements? (0=no=def, 1=yes)",
			   0, 0);
      stringstream ss;
      ss << "shapes-" << iTheta << "-" << nTheta << ".dat";
      parameters.addMember("outputFilename", &outputFilename, def,
			   "name of the outputFilename, if outputFlags > 0", ss.str().c_str());
  }

  //
  // set parameters
  //
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
    // check for valid iTheta (lower bound checked when generating Pset): 
    if (iTheta >= nTheta) {
      cerr << "iTheta value (" << iTheta << ") not compatible with nTheta (" << nTheta << ")" 
	   << endl;
      exit(1);
    }

    //
    // Open output file, if specified
    //
    std::ofstream ofs;
    if (outputFlags == 1) {
      ofs.open(outputFilename.c_str());
    }
    else if (outputFlags != 0) {
      cerr << "Invalid value for outputFlags:" << outputFlags << endl;
      exit(1);
    }

    cout << "# " << stringstuff::taggedCommandLine(argc, argv) << endl;
    parameters.dump(cout);
    if (outputFlags==1)  parameters.dump(ofs);

    //
    // First create the galaxy
    //
    SBProfile* unlensed=SBParse(sbGalaxy);
    // Measure half-light radius
    double ee50g = EnclosedFluxRadius(*unlensed);
    // Dilate if desired:
    if (rhalf > 0.) {
      cout << "# Raw galaxy half-light radius " << ee50g << endl;
      if (outputFlags==1)  ofs << "# Raw galaxy half-light radius " << ee50g << endl;
      double factor = rhalf / ee50g;
      Ellipse e(0.,0.,log(factor));
      SBProfile* sb2 = unlensed->distort(e);
      delete unlensed;
      unlensed = sb2;
      ee50g = rhalf;
    }
    cout << "# Galaxy EE50: " << ee50g << endl;
    if (outputFlags==1)  ofs << "# Galaxy EE50: " << ee50g << endl;

    //
    // Create PSF
    //
    SBProfile* psfraw=SBParse(sbPSF);
    double ee50psf = EnclosedFluxRadius(*psfraw);
    // Dilate if desired:
    if (rhalfPSF > 0.) {
      cout << "# Raw PSF half-light radius " << ee50psf << endl;
      if (outputFlags==1)  ofs << "# Raw PSF half-light radius " << ee50psf << endl;
      double factor = rhalfPSF / ee50psf;
      Ellipse e(0.,0.,log(factor));
      SBProfile* sb2 = psfraw->distort(e);
      delete psfraw;
      psfraw = sb2;
      ee50psf = rhalfPSF;
    }
    psfraw->setFlux(1.);	// Always unit normalized
    cout << "# PSF EE50: " << ee50psf << endl;
    if (outputFlags==1)  ofs << "# PSF EE50: " << ee50psf << endl;

    //
    // Draw a version of the object at zero rotation to decide noise
    // level that makes desired significance
    //
    double gflux=-1.;
    double imgRMS = -1.;
    double estSize = hypot( ee50g, ee50psf);
    Ellipse cleanBasis(0., 0., log(estSize));
    {
      SBConvolve obs(*unlensed, *psfraw);
      double ee30obs = EnclosedFluxRadius(obs, 0.3);
      double ee50obs = EnclosedFluxRadius(obs, 0.5);
      double ee70obs = EnclosedFluxRadius(obs, 0.7);
      cout << "# Observed EE30, 50, 70: "
	   << ee30obs << " " << ee50obs << " " << ee70obs << " " << endl;
      if (outputFlags==1)  ofs << "# Observed EE30, 50, 70: "
			       << ee30obs << " " << ee50obs << " " << ee70obs << " " << endl;
      gflux = obs.getFlux();
      cout << "# Galaxy flux: " << gflux << endl;
      if (outputFlags==1)  ofs << "# Galaxy flux: " << gflux << endl;
      // Set noise level from circular-aperture formula:
      imgRMS = 0.5*gflux / (sqrt(PI)*ee50obs*significance);
      cout << "# Pixel RMS: " << imgRMS << endl;
      if (outputFlags==1)  ofs << "# Pixel RMS: " << imgRMS << endl;
      
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
	//clean.shift(1,1);
	//FITSImage<>::writeToFITS("clean.fits", clean);
	exit(1);
      }
      cleanBasis = gl.getBasis();
      cleanBasis.setMu( cleanBasis.getMu() + log(dx));
      cout << "# Observed GL sigma: " << exp(cleanBasis.getMu()) << endl;
    if (outputFlags==1)  ofs << "# Observed GL sigma: " << exp(cleanBasis.getMu()) << endl;
      double f, varf;
      gl.b00(f, varf);
      double snr = (f / sqrt(varf)) * (significance / trialSN);
      cout << "# Observed GL S/N: " << snr << endl;
    if (outputFlags==1)  ofs << "# Observed GL S/N: " << snr << endl;

      // Calculate optimal detection significance
      double sigsum=0.;
      double fluxsum=0.;
      for (int iy=clean.YMin(); iy<=clean.YMax(); iy++)
	for (int ix=clean.XMin(); ix<=clean.XMax(); ix++) {
	  sigsum += clean(ix,iy)*clean(ix,iy);
	  fluxsum += clean(ix,iy);
	}
      cout << "# Observed ideal S/N: " << sqrt(sigsum)*dx/imgRMS << endl;
    if (outputFlags==1)  ofs << "# Observed ideal S/N: " << sqrt(sigsum)*dx/imgRMS << endl;

      // save the object as a FITS file
      //clean.shift(1,1);
      //FITSImage<>::writeToFITS("clean.fits", clean);
    }

    //
    // Make a PSFInfo structure - do a GL fit to get its basis
    //
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
      if (outputFlags==1)  ofs << "# PSF e and GL sigma: " << psfBasis.getS()
			       << " " << exp(psfBasis.getMu())
			       << endl;
    }

    //
    // begin: the iTheta "broken" loop
    //
    { 
      // things that we need local copies of, for parallel prossessing
      SBProfile* psfraw=SBParse(sbPSF);
      PSFInformation psfinfo(*psfraw, psfBasis);
      SBProfile* unlensed=SBParse(sbGalaxy);
      UnweightedShearEstimator seRawTheta;
      UnweightedShearEstimator seNoFixTheta;
      double sum_sn = 0.;

      UniformDeviate u;
      GaussianDeviate g(u);
      ostringstream oss;

      int nFailTheta = 0;	// Number of failures at this angle

      double theta = (2*PI * iTheta) / nTheta;
      
      oss << "# iTheta / nTheta = " << iTheta << " / " << nTheta << endl;
      oss << "# theta = " << theta / DEGREE << endl;

      SBProfile* gspin = unlensed->rotate(theta);
      SBProfile* gshear= gspin;
      if (dist1==0. && dist2==0.) 
	gspin=0;	// so we can know not to double-delete it
      else
	gshear = gspin->shear(dist1, dist2);

      Ellipse startBasis = cleanBasis;
      // Rotate the shear
      {
	double e1, e2;
	e1 = startBasis.getS().getE1();
	e2 = startBasis.getS().getE2();
	double e1r = cos(2*theta)*e1 -sin(2*theta) * e2;
	double e2r = sin(2*theta)*e1 +cos(2*theta) * e2;
	Shear srot(e1r, e2r);
	startBasis.setS(srot);
      }

      // Collect stats on this rotation:
      double sum1=0., sum2=0., sumsq1=0., sumsq2=0.;
      double sumFixEst = 0.;
      double sumSigma=0.;
      
      for (int iDither=0; iDither<dithersPerTheta; iDither++) {
	SBProfile* gshift = gshear->shift(u-0.5, u-0.5);
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
	      sci(ix,iy) += g*imgRMS;
	  FitExposure<> fe(sci, wt, xw, yw, 0);
	  FDNT<> fd(fe, psfinfo, startBasis, order);
	  fd.GLAll();
	  fd.setMaskSigma(6.);
	  bool success = fd.prepare();
	  Shear targetS = fd.getBasis().getS();  
	  double prob;
	  tmv::SymMatrix<double> covE(2,2);
	  if (success) {
	    try {  // catch FDNTError
	      targetS = fd.shape2(prob, covE);  // THIS IS THE ACTUAL MEASUREMENT!!
	    } catch (FDNTError &m) {
	      // Singularty flag has been set, success flag should catch it
	    }
	    // ??? Make flag mask an input parameter ???
	    success = !(fd.getFlags() & (DidNotConverge + Singularity
					 + OutOfBounds + TooLarge + UnderSampled
					 + TooElliptical + GLFailure));
	  }
	  if (!success) {
	    failNoise++;
	    nFailTheta++;
	    if ( ((realizationsPerDither > 100) && (failNoise > realizationsPerDither / 3)) 
		 || (failNoise > 4*realizationsPerDither)) {
	      // ?? Check limits above - quit if ~>50% of trials failing 
	      cerr << "Too many failures at one dither" << endl;
	      exit(1);
	    }
	    continue;
	  }

	  // Success: accumulate statistics here 
	  double f, varf;
	  fd.wtFlux(f, varf);
	  sum_sn += f / sqrt(varf);
	  sumSigma += exp(fd.getBasis().getMu());
	  iNoise++;
	  double egFix = fd.shrinkResponse(targetS);
	  sumFixEst += egFix;
	  seRawTheta.add(targetS, egFix);
	  seNoFixTheta.add(targetS, 0.);
	  double eta1, eta2;
	  targetS.getEta1Eta2(eta1,eta2);
	  sum1 += eta1;
	  sum2 += eta2;
	  sumsq1 += eta1*eta1;
	  sumsq2 += eta2*eta2;

	  // ...and save to output if necessary
	  if (outputFlags == 1) {
	    double ee30obs = EnclosedFluxRadius(sci, 0.3);
	    double ee50obs = EnclosedFluxRadius(sci, 0.5);
	    double ee70obs = EnclosedFluxRadius(sci, 0.7);
	    double g1, g2;
	    targetS.getG1G2(g1, g2);
	    ofs << g1 << " " << g2 << " " << eta1 << " " << eta2 << " "
		<< covE(0,0) << " " << covE(1,1) << " " << covE(0,1) << " "
		<< ee30obs << " " << ee50obs << " " << ee70obs << " "
		<< f << " " << sqrt(varf) << " " << exp(fd.getBasis().getMu()) << endl;
	  }

	} // end realization loop
      } // end dither loop

      if (gspin) delete gspin;
      delete gshear;

      // print out summary for this iTheta
      oss << "## trial counts info ##" << endl;
      int N = dithersPerTheta * realizationsPerDither;
      oss << "# iTheta, theta(deg), N, nFailTheta, avg(S/N)" << endl;
      oss << iTheta << " " << theta/DEGREE << " " << N << " " << nFailTheta << " " << sum_sn/N
	  << endl;

      oss << "## seMeans info ##" << endl;
      oss << "# sumeta1, sumeta2, sumetasq1, sumetasq2, N, sumFixEst, theta(deg)" << endl;
      oss << sum1 << " " << sum2 << " " << sumsq1 << " " << sumsq2 << " " 
	  << N << " " << sumFixEst << " " << theta/DEGREE << endl;

      oss << "## seRaw info ##" << endl;
      oss << seRawTheta.getLine();

      oss << "## seNoFix info ##" << endl;
      oss << seNoFixTheta.getLine();

      cout << oss.str();
      if (outputFlags==1)  ofs << oss.str();
    } // end iTheta "loop"

  } catch (tmv::Error& m) {
    cerr << m << endl;
  } catch (std::runtime_error &m) {
    quit(m,1);
  }
}



      /*/ Print out statistics for this theta
      double N = dithersPerTheta * realizationsPerDither;
      double mean1 = sum1 / N;
      double mean2 = sum2 / N;
      double sig1 = sqrt(sumsq1/N-mean1*mean1);
      double sig2 = sqrt(sumsq2/N-mean2*mean2);
      double egFix = sumFixEst / N;
      oss << setw(2) << iTheta
	  << " " << setw(8) << mean1
	  << " " << setw(8) << mean2
	  << " +- " << setw(8) << sig1
	  << " " << setw(8) << sig2
	  << " " << (1.*nFailTheta) / (nFailTheta
				       +dithersPerTheta*realizationsPerDither)
	  << " " << sumSigma / N
	  << " " << egFix
	  << endl;
      double e1derotate = mean1 * cos(2*theta) + mean2 * sin(2*theta);
      double e2derotate = -mean1 * sin(2*theta) + mean2 * cos(2*theta);
      /*/

      /*/ Add into overall estimator:

    // Initialize statistics: one that takes mean of each ring pos,
    // one that takes all galaxies together
    double sumFixAll = 0.;
    double sumEta1der = 0.;
    double sumEta2der = 0.;


      {
	Shear S;
	S.setEta1Eta2(mean1, mean2);
	seMeans.add(S, egFix, sig1/sqrt(N-1), sig2/sqrt(N-1));

	seRaw.addEnsemble(seRawTheta);
	seNoFix.addEnsemble(seNoFixTheta);

	sumEta1der += e1derotate;
	sumEta2der += e2derotate;
	sumFixAll += egFix;
	
      }
      /*/

    /*/
    Shear S(seMeans);
    double e1, e2, sig1, sig2;
    S.getEta1Eta2(e1,e2);
    seMeans.sigmaE(sig1,sig2, false);

    cout << fixed << setprecision(6);
    cout << "Means: " << e1 << " +- " << sig1
	 << " " << e2 << " +- " << sig2 
	 << setprecision(5)
	 << " derotate " << sumEta1der/nTheta
	 << " " << sumEta2der / nTheta
	 << setprecision(4)
	 << " egfix " << sumFixAll / nTheta
	 << endl;

    cout << fixed << setprecision(6);
    cout << "Raw: " << (Shear) seRaw
	 << endl;
    cout << "NoFix: " << (Shear) seNoFix
	 << endl;
    /*/
