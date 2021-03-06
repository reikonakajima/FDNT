// $Id: EnclosedFluxRadius.h,v 1.2 2011/07/20 15:44:35 dgru Exp $
// Report the radius enclosing desired flux fraction
// for an Image or an SBProfile.
#ifndef ENCLOSEDFLUXRADIUS_H
#define ENCLOSEDFLUXRADIUS_H
#include "Image.h"
#include "SBProfile.h"

double EnclosedFluxRadius(Image<> img, double enclosedFraction=0.5);
double EnclosedFluxRadius(Image<> img, double xc, double yc, double enclosedFractionFlux);

// Find EEradius for an SBProfile
double EnclosedFluxRadius(const sbp::SBProfile& sbp, double enclosedFraction=0.5);

#endif // ENCLOSEDFLUXRADIUS_H
