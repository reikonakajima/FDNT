//
// === ColorGrad.h ===
//
// 20121212 RN  : initial draft
//
// 
// ColorGradGaussian class:
// 
// Reads in E. Semboloni's color gradient PSF data file
// [format:
//   (1) wavelength  (2) gaussian width sigma  (3) normalization (SED of star)
// ]
// and generate a SBProfile (actually of an SBAdd class) describing the PSF.
//
#include <vector>
#include <list>
#include <stdio.h>
#include <stdlib.h>
#include "Std.h"
#include "StringStuff.h"
#include "SBProfile.h"
#include "FITSImage.h"

#ifndef COLORGRAD_H
#define COLORGRAD_H

using namespace std;
using namespace sbp;

class ColorGradGaussian {
 public:
  ColorGradGaussian(string data_file_name);
  SBAdd& getSBProfile() {return psf;}
 private:
  SBAdd psf;
};

#endif    // COLORGRAD_H
