#ifndef RUNFDNT_H
#define RUNFDNT_H

#include "FDNT.h"
#include "StringStuff.h"
#include "FITSImage.h"
#include <fstream>
#include "ShearEstimator.h"
#include "EnclosedFluxRadius.h"
#include "Image.h"

/// future projects, to include PSF model and astrometry
//#include "PSFEx.h"
#include "Astrometry.h"
//#include "SCAMPMap.h"


namespace galsim {
namespace fdnt{

// a simple test function for the module.  returns the string 'hello, world'
char const* greet();


/*
 * RunFDNT: the c++ function that interfaces the Python layer.
 */
template <typename T>
int RunFDNT(const Image<T>& gal_image, const Image<T>& psf_image, const Image<T>& weight_image,
	    double x_pix, double y_pix,
	    double a_wc, double b_wc, double pa_wc,
	    double r_pix, // FLUX_RADIUS in pixels, not wcs
	    double ee50psf,
	    double bg,
	    int order, double sky);


}
}


#endif // RUNFDNT_H
