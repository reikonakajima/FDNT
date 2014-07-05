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


namespace fdnt{

  // a simple test function for the module.  returns the string 'hello, world'
  char const* greet();


  /*
   * @brief C++ Struct containing information about the shape of the object.
   *
   * This representation contains information about (i) observed shapes, (ii) PSF shapes,
   * and (iii) shape estimates after PSF correction, via the FDNT method.
   *
   * Observed shapes (convolved with PSF) have |e|< 1 or |g|< 1, while the PSF corrected shapes
   * may exceed 1 in its absolute value.
   *
   */
  struct FDNTShapeData {
    /// @brief galsim::Bounds object describing the image of the object
    Bounds<int> image_bounds;

    //
    // (i) Observed galaxy information:
    //

    /// @brief Status after measuring adaptive moments (via GLSimple).  -1 if not set.
    int observed_flags;

    /// @brief Distortion component e1 of the observed shape.  -10. if not measured.
    double observed_e1;

    /// @brief Distortion component e2 of the observed shape.  -10. if not measured.
    double observed_e2;

    /// @brief Size sigma=exp(mu) of observed galaxy, currently (!) in units of pixels
    float observed_sigma;

    /// @brief Total image intensity of the best-fit elliptical Gaussian from adaptive moments;
    /// note that this is not flux, since flux = (total image intensity)*(pixel scale)^2
    double observed_b00;

    /// @brief Variance in the total image intensity of the best-fit elliptical Gaussian.
    double observed_b00_var;

    /// @brief The weighted radial fourth moment of the image.  -1. if not measured.
    double observed_b22;

    /// @brief Centroid of best-fit elliptical Gaussian
    Position<double> observed_centroid;

    //
    // (ii) Estmiated PSF information:
    //

    /// @brief Status after measuring PSF adaptive moments (via GLSimple).  -1 if not set.
    int psf_flags;

    /// @brief Distortion component e1, representing the PSF shape.  -10. if not measured.
    double psf_e1;

    /// @brief Distortion component e2, representing the PSF shape.  -10. if not measured.
    double psf_e2;

    /// @brief Size sigma=exp(mu) of the PSF, currently in units of pixels
    float psf_sigma;

    /// @brief The GL polynomial order to which the PSF was fit to.  -1 if not measured.
    int psf_order;

    /// @brief Total image intensity of the best-fit elliptical Gaussian from adaptive moments;
    /// note that we operate in the pixel scale.  -1. if not measured.
    double psf_b00;

    /// @brief Variance in the psf_amp.  -1. if not measured.
    double psf_b00_var;

    /// @brief The weighted radial fourth moment of the image.  -1. if not measured.
    double psf_b22;

    /// @brief Chi squared value of the PSF fit.  0 if not measured.
    double psf_chisq;

    /// @brief Degree of freedom in the PSF fit.  0 if not measured.
    int psf_DOF;

    //
    // (iii) PSF-corrected galaxy information:
    //

    /// @brief Status after PSF-corrected shape measurements (via FDNT).  -1 if not set.
    int intrinsic_flags;

    /// @brief Estimated e1 after correcting for effects of the PSF.  -10 if not measured.
    float intrinsic_e1;

    /// @brief Estimated e2 after correcting for effects of the PSF.  -10 if not measured.
    float intrinsic_e2;

    /// @brief Shape measurement variance in e1
    float intrinsic_e1_var;

    /// @brief Shape measurement variance in e2
    float intrinsic_e2_var;

    /// @brief Shape measurement covariance in e1 and e2
    float intrinsic_e1e2_covar;

    /// @brief Size sigma=exp(mu) of the intrinsic galaxy, currently in units of pixels.
    float intrinsic_sigma;

    /// @brief Bias in the responsivity due to ellipticity gradient.  -10. if not set.
    double shrink_response;

    /// @brief Evaluation count for the FDNT measurement.
    int evaluation_count;

    /// @brief Trial ellipticity count for the FDNT measurement.
    int e_trial_count;

    /// @brief Resolution factor R_2; 0 indicates object is consistent with a PSF, 1 indicates
    /// perfect resolution; default -1
    float resolution_factor;

    /// @brief A string containing any error messages from the attempted measurements, to
    /// facilitate proper error handling in both C++ and python
    std::string error_message;

    /// @brief Constructor, setting defaults
    FDNTShapeData() : image_bounds(Bounds<int>()), observed_flags(-1), observed_e1(-10.),
      observed_e2(-10.), observed_sigma(-1.), observed_b00(-1.), observed_b00_var(-1.),
      observed_b22(-1.), observed_centroid(Position<double>(0.,0.)),

      psf_flags(-1), psf_e1(-10.), psf_e2(-10.), psf_sigma(-1.), psf_order(-1),
      psf_b00(-1.), psf_b00_var(-1.), psf_b22(-1.), psf_chisq(0.), psf_DOF(0),

      intrinsic_flags(-1), intrinsic_e1(-10.), intrinsic_e2(-10.), intrinsic_e1_var(0.),
      intrinsic_e2_var(0.), intrinsic_e1e2_covar(0.), intrinsic_sigma(-1.),
      shrink_response(-10.), evaluation_count(0), e_trial_count(0),
      resolution_factor(-1.), error_message("")
      {}

  };



  /*
   * @brief The C++ shape measurement function that interfaces the Python layer.
   * The input parameter takes an initial guess values from, e.g., SExtractor.
   *
   * Note that currently the WCS information from thie Image<> class is not used, which means
   * that all initial guess information provided (x_wc, y_wc, a_wc, b_wc, pa_wc) must be
   * in pixel units.  When the WCS is implemented, the x_wc and y_wc should change to
   * x_wcs and y_wcs.  sigma_pix will remain in pixel units, regardless of whether the WCS is
   * implemented or not.
   *
   * The FDNT will run on a postage stamp that must be (minimumStampSigma * sigma_pix) in
   * dimention or larger.  FDNT will throw if the image size is smaller than this.
   */
  template <typename T>
    FDNTShapeData RunFDNT(const Image<T>& gal_image, const Image<T>& psf_image,
			   const Image<T>& weight_image,
			   double x_wc, double y_wc,
			   double a_wc, double b_wc, double pa_wc,
			   double sigma_pix, // FLUX_RADIUS in pixels, not wcs
			   double ee50psf,
			   double bg,
			   int order, double sky);

}


#endif // RUNFDNT_H
