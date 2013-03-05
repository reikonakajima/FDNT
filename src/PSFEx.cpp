// $Id: PSFEx.cpp,v 1.10 2012/08/13 17:51:29 dgru Exp $
#include "PSFEx.h"
#include "Interpolant.h"
#include <fstream>

using namespace sbp;

// Create the interpolant assumed to be used by PSFEx:
// Lanzcos3, flux-conserving, tolerance=1e-3
const Lanczos lanczos3(3, true, 1e-3);
InterpolantXY 
PSFExModel::lanczos3_2d(lanczos3);
// The cubic and quintic are the interpolants to use in k space:
const Cubic cubic(1e-4);
const Quintic quintic(1e-4);
InterpolantXY 
PSFExModel::cubic_2d(cubic);
InterpolantXY 
PSFExModel::quintic_2d(quintic);

void
PSFExModel::fieldPosition(double x, double y, int maxorder) {
// use only up to maxorder polynomial order from the model (unless maxorder is <= 0, then use all available monomials)

  if(xOrder.size()<2) // catch constant PSF
  {
   DVector wts(1);
   wts[0]=1.;
   sbp->setWeights(wts);
   return;
  }
  double xscale = (x-polzero[0])/polscal[0];
  double yscale = (y-polzero[1])/polscal[1];
  DVector wts(xOrder.size());
  for (int i=0; i<xOrder.size(); i++) {
    if(maxorder>=0 && xOrder[i]+yOrder[i]>maxorder)
      wts[i]=0.;
    else
      wts[i] = pow(xscale, xOrder[i]) * pow(yscale, yOrder[i]);
  }
  sbp->setWeights(wts);
}

PSFExModel::PSFExModel(const char *filename, PSFExFormat format): sbp(0), fwhm(0.), dx(0.),
								  Nx(0), polnaxis(0)
{
  if(format!=TEXT)
    throw PSFExError("PSFEx only implemented for TEXT input format");

  ifstream in(filename);
	
  if(!in.good())
    FormatAndThrow<PSFExError>() << "Could not open PSFEx model file " << filename;
	
  string buf;
  bool headerread=false; // flag set true once header has been read completely
  bool polnaxisread=false;
  bool polngrpread=false;
  bool psfnaxisread=false;
  int polngrp = 0;
  int polnaxis = 0;
  int psfnaxis = 0;
  double chi2;	// discarded after reading.
  vector<int> psfaxis;
  
  do {
    in >> buf;
    if (!buf.compare("LOADED"))	{
      int junk; in >> junk;
      // Not used
    } else if(!buf.compare("ACCEPTED")) {
      int junk; in >> junk;
      // Not used
    } else if(!buf.compare("POLNAXIS"))	{
      if( !(in >> polnaxis) || (polnaxis!=2 && polnaxis) )
	FormatAndThrow<PSFExError> () << "PSFEx only accepts POLNAXIS==0 or 2, read " << polnaxis;
      polgrp.resize(polnaxis, 0);
      polname.resize(polnaxis, "UNKNOWN");
      polzero.resize(polnaxis, 0.);
      polscal.resize(polnaxis, 1.);
      polnaxisread=true;

    } else if(!buf.compare(0,6,"POLGRP")) {
      int i = atoi(buf.substr(6).c_str())-1;
      if (!polnaxisread || i<0 || i>=polnaxis || !(in >> polgrp[i]))
	  throw PSFExError("PSFEx read bad or out-of-order POLGRP");
      // Note not using polgrp right now.

    } else if(!buf.compare(0,7,"POLNAME")) {
      int i = atoi(buf.substr(7).c_str())-1;
      if (!polnaxisread || i<0 || i>=polnaxis || !(in >> polname[i]))
	  throw PSFExError("PSFEx read bad or out-of-order POLNAME");
      // Note polname is not used at present

    } else if(!buf.compare(0,7,"POLZERO")) {
      int i = atoi(buf.substr(7).c_str())-1;
      if (!polnaxisread || i<0 || i>=polnaxis || !(in >> polzero[i]))
	throw PSFExError("PSFEx read bad or out-of-order POLZERO");

    } else if(!buf.compare(0,7,"POLSCAL")) {
      int i = atoi(buf.substr(7).c_str())-1;
      if (!polnaxisread || i<0 || i>=polnaxis || !(in >> polscal[i]))
	throw PSFExError("PSFEx read bad or out-of-order POLSCAL");

    } else if(!buf.compare("POLNGRP")) {
      if (!(in >> polngrp))
	throw PSFExError("PSFEx could not read POLNGRP");
      if (polngrp>1) 
	FormatAndThrow<PSFExError> () << "PSFEx can only use POLNGRP<=1, read " << polngrp;

      poldeg.resize(polngrp, 0);
      polngrpread = true;

    } else if(!buf.compare(0,6,"POLDEG")) {
      int i = atoi(buf.substr(6).c_str())-1;
      if (!polngrpread || i<0 || i>=polngrp || !(in >> poldeg[i]) )
	throw PSFExError("PSFEx read bad or out-of-order POLDEG");

    } else if (!buf.compare("PSF_FWHM")) {
      if (!(in >> fwhm))
	throw PSFExError("PSFEx could not read PSF_FWHM");

    } else if (!buf.compare("PSF_SAMP")) {
      if (!(in >> dx))
	throw PSFExError("PSFEx could not read PSF_SAMP");

    } else if (!buf.compare("CHI2")) {
      if (!(in >> chi2))
	throw PSFExError("PSFEx could not read CHI2");
      // Note this is not retained.

    }  else if (!buf.compare("PSFNAXIS")) {
      if (!(in >> psfnaxis)) 
	throw PSFExError("PSFEx could not read PSFNAXIS");
      if (psfnaxis!=3) 
	FormatAndThrow<PSFExError> () << "PSFEx can only use PSFNAXIS==3, read " << psfnaxis;

      psfaxis.resize(psfnaxis, 0);
      psfnaxisread=true;

    } else if (!buf.compare(0,7,"PSFAXIS")) {
      int i = atoi(buf.substr(7).c_str())-1;
      if(!psfnaxisread || i < 0 || i>=psfnaxis || !(in >> psfaxis[i]) || psfaxis[i]<1)
	throw PSFExError("PSFEx read bad or out-of-order PSFAXIS");

    } else if (!buf.compare("DATA")) {
      headerread=true;
      // DATA keyword marks end of header information

    } else {
      FormatAndThrow<PSFExError> () << "PSFEx read unknown header word " << buf;
    }

  } while(!headerread && !in.eof());
  
  if ( !headerread || !psfnaxisread || !polngrpread || !polnaxisread 
     || psfnaxis!=3 || psfaxis[0]<=0 || psfaxis[1]<=0 || psfaxis[2]<=0)
    FormatAndThrow<PSFExError>() << "PSFEx is missing header information in " << filename;

  // Make sure that psfaxis[2] matches the number of polynomial terms needed:
  xOrder.clear();
  yOrder.clear();
  if(polngrp>0)
   order = poldeg[0];
  else
   order = 0;
  for(int ny=0; ny<=order; ny++) { 
    for (int nx = 0; nx<=order-ny; nx++) {
      xOrder.push_back(nx);
      yOrder.push_back(ny);
    }
  }
  if (psfaxis[2] != xOrder.size())
    FormatAndThrow<PSFExError> () << "PSFEx polynomial term counts do not match: psfaxis[2]="
				  << psfaxis[2]
				  << " but order " << order
				  << " requires " << xOrder.size();

  // Create our SBPixel now:
  Nx = MAX(psfaxis[0], psfaxis[1]);
  sbp = new SBPixel(Nx, dx, lanczos3_2d, xOrder.size());
  // Now read the kernel elements
  for(int k=0; k<xOrder.size(); k++) {
    for(int j=0; j<psfaxis[1]; j++) {
      for(int i=0; i<psfaxis[0]; i++) {
	double value;
	if (!(in >> value)) 
	  throw PSFExError("PSFEx bad or missing data value");
	// Put values into the SBPixel arrays.  Assume that center ( 0,0 in SBPixel)
	// is at the (N/2+1)th element of the input array, which is N/2 for zero-indexed
	// counters
	sbp->setPixel(value, i-psfaxis[0]/2, j-psfaxis[1]/2, k);
      }
    }
  }
}

bool PSFExModel::selfTest(double tolerance, int grid, int order) {
// very simple self test of the PSF model, where we demand that the pixel with the least flux has more than -tolerance*brightest_pixel flux at some positions sampled over the field
  int n = sbp->getNin();
  
  for (int ix=0; ix<grid; ix++)
    for (int iy=0; iy<grid; iy++)
    {
      if(xOrder.size()>1)
       fieldPosition(polzero[0]+(double(ix)/(grid-1.)-0.5)*polscal[0],polzero[1]+(double(iy)/(grid-1.)-0.5)*polscal[1], order);
      else
       fieldPosition(0,0,order);
      double pmin=0.;
      double pmax=0.;
      double pminx=0.;
      double pminy=0.;
      double pmaxx=0.;
      double pmaxy=0.;
      for(int py=-n/2; py<n/2; py++)
	for(int px=-n/2; px<n/2; px++)
	{
	  double p = sbp->getPixel(px,py,-1);
	  if(p>pmax) { pmax=p; pmaxx=px; pmaxy=py; }
	  if(p<pmin) { pmin=p; pminx=px; pminy=py; }
        }
      if(-pmin>tolerance*pmax)
      {
	cerr << "PSFEx self-test failed at pmax=" << pmax << " in pixel " << pmaxx << "," << pmaxy << "; pmin=" << pmin << " in pixel " << pminx << "," << pminy << endl;
	cerr << "near field position " << polzero[0]+(double(ix)/(grid-1.)-0.5)*polscal[0] << "," << polzero[1]+(double(iy)/(grid-1.)-0.5)*polscal[1] << endl;
	cerr << sbp->getWeights();
	return false;
      }
    }
   return true;
}

PSFExModel::~PSFExModel() {
  if (sbp) delete sbp;
}

