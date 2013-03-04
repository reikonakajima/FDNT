// $Id: EnclosedFluxRadius.h,v 1.2 2011/07/20 15:44:35 dgru Exp $
// Report the radius enclosing desired flux fraction
// for an Image or an SBProfile.

#include "Image.h"
#include "SBProfile.h"

extern double EnclosedFluxRadius(Image<> img, double enclosedFraction=0.5);

// Find EEradius for an SBProfile
extern double EnclosedFluxRadius(const sbp::SBProfile& sbp, double enclosedFraction=0.5);

  
