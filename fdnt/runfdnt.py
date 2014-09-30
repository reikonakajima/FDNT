""" @file runfdnt.py
Routines for shape measurements with the FDNT method.

This file contains the python interface to C++ routines for shape measurement, given the galaxy,
PSF and mask images (reference: Bernstein 2010, Bernstein & Jarvis 2002).

"""


from . import _fdnt
import fdntimage
import bounds
import galsim
import numpy as np

class FDNTShapeData(object):
    """A class to contain the outputs of using the FDNT shape and moments measurement routines.

    At the C++ level, there is a container for the outputs of the FDNT shape measurement routines.
    The FDNTShapeData class is the analogous object at the python level.  It contains the following
    information about moment measurement (from either FDNTEstimateShear() or GLMoments()):

    - image_bounds: a BoundsI object describing the image.

    - observed_flags: the status flag resulting from GLMoments measurement; -1 indicates no
      attempt to measure.

    - observed_e1: distortion component e1 of the observed shape based on Gauss-Laguerre fits.

    - observed_e2: distortion component e2 of the observed shape based on Gauss-Laguerre fits.

    - observed_e1_var: variance in distortion component e1

    - observed_e2_var: variance in distortion component e2

    - observed_e1e2_covar: covariance in distortion components e1 and e2

    - observed_sigma: size sigma of the best-fit Gaussian, in units of pixels; -1 if
      not measured.

    - observed_b00: total image intensity for best-fit elliptical Gaussian from the best GLMoment
      fit.  Normally, this field is simply equal to the image flux (for objects that follow a
      Gaussian light distribution, otherwise it is something approximating the flux).  However,
      if the image was drawn using `drawImage(method='sb')` then observed_amp relates to the flux
      via flux = (observed_amp)*(pixel scale)^2.

    - observed_b00_var: variance in the observed_amp value.

    - observed_b22: the weighted radial fourth moment of the image.

    - observed_centroid: a PositionD object representing the centroid based on GLMoment fit.
      The convention for centroids is such that the center of the lower-left pixel is
      (0,0).

    - observed_significance: signal-to-noise, defined as the ratio
      [observed_b00 / sqrt(observed_b00_var)].

    If FDNTEstimateShear() was used, then the following fields related to PSF-corrected shape will
    also be populated:

    - psf_flags: Status after measuring PSF via GLMoments().  -1 if not set.

    - psf_e1, psf_e2: Distortion components e1/e2 of the PSF shape via GLMoments().
      -10. if not measured.

    - psf_sigma: Size sigma=exp(mu) of the PSF, currently in units of pixels.

    - psf_order: The GL polynomial order to which the PSF was fit to.  -1 if not measured.

    - psf_b00: Total image intensity of the best-fit elliptical Gaussian from GLMoments().
      Units [currently] in the pixel scale.  -1. if not measured.

    - psf_b00_var: Variance in psf_b00.  -1. if not measured.

    - psf_b22: The weighted radial fourth moment of the image.  -1. if not measured.

    - psf_chisq: Chi squared value of the GLMoment fit of the PSF.  0 if not measured.

    - psf_DOF: Degree of freedom in the GLMoment fit of the PSF.  0 if not measured.

    - intrinsic_flags: the status flag resulting from FDNT measurements; -1 indicates no attempt to
      measure, 0 indicates success.

    - instrinsic_e1, instrinsic_e2: doubles representing the estimated shear after removing
      the effects of the PSF.  e1/e2 will differ from the default values of -10 if a meaurement
      is made.  e1/e2 corresponds to distortion defined by e = (a^2-b^2)/(a^2+b^2), not shear
      g = (a-b)/(a+b), where a and b are the major and minor axis, respectively.

    - instrinsic_e1_var, instrinsic_e2_var, intrinsic_e1e2_covar: doubles representing the
      estimated covariance elements of the shear after removing the effects of the PSF.
      -10. if a meaurement was not made.

    - intrinsic_sigma: Size sigma=exp(mu) of the intrinsic galaxy, currently in units of pixels.

    - shrink_response: Bias in the responsivity due to ellipticity gradient.  -10. if not set.

    - evaluation_count: Evaluation count for the FDNT measurement.

    - e_trial_count: Trial ellipticity count for the FDNT measurement.

    - resolution_factor: Resolution factor R_2; 0 indicates object is consistent with a PSF,
      1 indicates perfect resolution; default -1

    - error_message: A string containing any error messages from the attempted measurements, to
      facilitate proper error handling in both C++ and python


    The FDNTShapeData object can be initialized completely empty, or can be returned from the
    routines that measure object moments (FindAdaptiveMom()) and carry out PSF correction
    (EstimateShear()).
    """
    def __init__(self, *args):
        # arg checking: require either a _FDNTShapeData, or nothing
        if len(args) > 1:
            raise TypeError("Too many arguments to initialize FDNTShapeData!")
        elif len(args) == 1:
            if not isinstance(args[0], _fdnt._FDNTShapeData):
                raise TypeError("Argument to initialize FDNTShapeData must be a _FDNTShapeData!")
            self.image_bounds = args[0].image_bounds
            self.observed_flags = args[0].observed_flags
            self.observed_e1 = args[0].observed_e1
            self.observed_e2 = args[0].observed_e2
            self.observed_e1_var = args[0].observed_e1_var
            self.observed_e2_var = args[0].observed_e2_var
            self.observed_e1e2_covar = args[0].observed_e1e2_covar
            self.observed_sigma = args[0].observed_sigma
            self.observed_b00 = args[0].observed_b00
            self.observed_b00_var = args[0].observed_b00_var
            self.observed_b22 = args[0].observed_b22
            self.observed_centroid = args[0].observed_centroid
            self.observed_significance = args[0].observedSignificance()

            self.psf_flags = args[0].psf_flags
            self.psf_e1 = args[0].psf_e1
            self.psf_e2 = args[0].psf_e2
            self.psf_sigma = args[0].psf_sigma
            self.psf_order = args[0].psf_order
            self.psf_b00 = args[0].psf_b00
            self.psf_b00_var = args[0].psf_b00_var
            self.psf_b22 = args[0].psf_b22
            self.psf_chisq = args[0].psf_chisq
            self.psf_DOF = args[0].psf_DOF

            self.intrinsic_flags = args[0].intrinsic_flags
            self.intrinsic_e1 = args[0].intrinsic_e1
            self.intrinsic_e2 = args[0].intrinsic_e2
            self.intrinsic_e1_var = args[0].intrinsic_e1_var
            self.intrinsic_e2_var = args[0].intrinsic_e2_var
            self.intrinsic_e1e2_covar = args[0].intrinsic_e1e2_covar
            self.intrinsic_sigma = args[0].intrinsic_sigma
            self.shrink_response = args[0].shrink_response
            self.evaluation_count = args[0].evaluation_count
            self.e_trial_count = args[0].e_trial_count
            self.resolution_factor = args[0].resolution_factor
            self.error_message = args[0].error_message

        else:
            self.image_bounds = _fdnt.BoundsI()
            self.observed_flags = -1
            self.observed_e1 = -10.
            self.observed_e2 = -10.
            self.observed_e1_var = -1.
            self.observed_e2_var = -1.
            self.observed_e1e2_covar = -1.
            self.observed_sigma = -1.0
            self.observed_b00 = -1.0
            self.observed_b00_var = -1.0
            self.observed_rho4 = -1.0
            self.observed_centroid = _fdnt.PositionD()
            self.observed_significance = 0.

            self.psf_flags = -1
            self.psf_e1 = -10.
            self.psf_e2 = -10.
            self.psf_sigma = -1.
            self.psf_order = -1
            self.psf_b00 = args[0].psf_b00
            self.psf_b00_var = args[0].psf_b00_var
            self.psf_b22 = args[0].psf_b22
            self.psf_chisq = args[0].psf_chisq
            self.psf_DOF = args[0].psf_DOF

            self.intrinsic_flags = -1
            self.intrinsic_e1 = -10.
            self.intrinsic_e2 = -10.
            self.intrinsic_e1_var = 0.
            self.intrinsic_e2_var = 0.
            self.intrinsic_e1e2_covar = 0.
            self.intrinsic_sigma = -1.
            self.shrink_response = -10.
            self.evaluation_count = 0
            self.e_trial_count = 0
            self.resolution_factor = -1.0
            self.error_message = ""


# A helper function for taking input weight and badpix Images, and returning a weight Image in the
# format that the C++ functions want
def _convertMask(image, weight=None, badpix=None):
    """Convert from input weight and badpix images to a single mask image needed by C++ functions.

    This is used by RunFDNT().
    """
    # if no weight image was supplied, make a float array (same size as gal image) filled with 1's
    if weight == None:

        # convert types from galsim.BoundsI to fdnt.BoundsI
        b = _fdnt.BoundsI(image.bounds.xmin, image.bounds.xmax, image.bounds.ymin, image.bounds.ymax)
        # create mask image
        mask = _fdnt.FDNTImageF(b, 1.)

    else:
        # if weight image was supplied, check if it has the right bounds and is non-negative
        if weight.bounds != image.bounds:
            raise ValueError("Weight image does not have same bounds as the input Image!")

        # also make sure there are no negative values
        import numpy as np
        if np.any(weight.array < 0) == True:
            raise ValueError("Weight image cannot contain negative values!")

        # if weight is an ImageI, then we can use it as the mask image:
        if weight.dtype == np.int32:
            if not badpix:
                mask = weight
            else:
                # If we need to mask bad pixels, we'll need a copy anyway.
                mask = _fdnt.FDNTImageF(weight)

        # otherwise, we need to convert it to the right type
        else:
            mask = galsim.ImageI(bounds=image.bounds, init_value=0)
            mask.array[weight.array > 0.] = 1

    # if badpix image was supplied, identify the nonzero (bad) pixels and set them to zero in weight
    # image; also check bounds
    if badpix != None:
        if badpix.bounds != image.bounds:
            raise ValueError("Badpix image does not have the same bounds as the input Image!")
        import numpy as np
        mask.array[badpix.array != 0] = 0

    # if no pixels are used, raise an exception
    if mask.array.sum() == 0:
        raise RuntimeError("No pixels are being used!")

    # finally, return the Image for the weight map
    return mask   ## all this copying is dangerous---if one is modified, so will the other!


def RunFDNT(gal_image, PSF_image, guess_x_wc, guess_y_wc,
            guess_sig_gal_pix, guess_sig_PSF_pix, guess_a_wc, guess_b_wc, guess_pa_wc,
            weight=None, order=0, bg=0., sky=0., badpix=None):
    """Carry out Fourier Domain Null Test PSF-corrected shape measurement routines.

    Example usage
    -------------

    Typical application to a single object:

        >>> galaxy = galsim.Gaussian(flux = 1.0, sigma = 1.0)
        >>> galaxy = galaxy.shear(g1=0.05, g2=0.0)  # shears the Gaussian by (0.05, 0) using the
        >>>                                         # |g| = (a - b)/(a + b) definition
        >>> psf = galsim.Kolmogorov(flux = 1.0, fwhm = 0.7)
        >>> final = galsim.Convolve([galaxy, psf])
        >>> final_image = final.drawImage(dx = 0.2)
        >>> final_epsf_image = psf.drawImage(dx = 0.2)
        >>> result = fdnt.RunFDNT(final_image, final_epsf_image)

    After running the above code, `result.observed_shape` ["shape" = distortion, the 
    (a^2 - b^2)/(a^2 + b^2) definition of ellipticity] is
    `(0.XX, 0.XX)` and `result.corrected_e1`, `result_corrected_e2` are
    `(0.XX, 0.XX)`, compared with the  expected `(0.09975, 0)` for a perfect
    PSF correction method.

    @param gal_image       The Image of the galaxy being measured.
    @param PSF_image       The Image for the PSF.
    @param weight          The optional weight image for the galaxy being measured.  Can be an int
                           or a float array.  Currently, GalSim does not account for the variation
                           in non-zero weights, i.e., a weight map is converted to an image with 0
                           and 1 for pixels that are not and are used.  Full use of spatial
                           variation in non-zero weights will be included in a future version of
                           the code.
    @param guess_x_wc      An initial guess for the x component of the object centroid (useful in
                           case it is not located at the center, which is the default
                           assumption).  The convention for centroids is such that the center of
                           the lower-left pixel is (0,0). [default: gal_image.trueCenter().x]
    @param guess_y_wc      An initial guess for the y component of the object centroid (useful in
                           case it is not located at the center, which is the default
                           assumption).  The convention for centroids is such that the center of
                           the lower-left pixel is (0,0). [default: gal_image.trueCenter().y]
    @param guess_sig_gal   Optional argument with an initial guess for the Gaussian sigma of the
                           galaxy (in pixels).
    @param guess_sig_PSF   Optional argument with an initial guess for the Gaussian sigma of the
                           PSF (in pixels).
    @param order           The optional GL-fit order, when GL fits are used to fill masked-out
                           pixels.  [default: 0]
    @param bg              The optional background level. [default: 0.]
    @param sky             The optional sky level.  [default: 0.]

    @returns a ShapeData object containing the results of shape measurement.
    """
    # prepare inputs to C++ routines:
    # _RunFDNT(img::Image<float> gal_image, img::Image<float> psf_image,
    #          img::Image<float> weight_image,
    #          double x_pix, double y_pix, double a_wc, double b_wc, double pa_wc,
    #          double r_pix, double ee50psf, double bg, int order, double sky)
    ## convert image formats
    gal_fdnt_image = _fdnt.FDNTImageF(gal_image.array, gal_image.xmin, gal_image.ymin)
    PSF_fdnt_image = _fdnt.FDNTImageF(PSF_image.array, PSF_image.xmin, PSF_image.ymin)
    weight_fdnt_image = _convertMask(gal_image, weight=weight, badpix=badpix)
    ## convert int to float
    guess_x_wc = float(guess_x_wc)   # centroids often specified by integers (pixels)
    guess_y_wc = float(guess_y_wc)

    try:
        result = _fdnt._RunFDNT(gal_fdnt_image, PSF_fdnt_image, weight_fdnt_image,
                                x_pix=guess_x_wc, y_pix=guess_y_wc,  # currently wc == pix
                                a_wc=guess_a_wc, b_wc=guess_b_wc, pa_wc=guess_pa_wc,
                                r_pix=guess_sig_gal_pix, ee50psf=guess_sig_PSF_pix,
                                bg=bg, order=order, sky=sky)

    except RuntimeError as err:
        raise RuntimeError

    return FDNTShapeData(result)



def GLMoments(gal_image, guess_x_wc, guess_y_wc,
              guess_sig_gal_pix, guess_a_b_pa=None, guess_g1g2=None,
              weight=None, order=0, bg=0., sky=0., badpix=None):
    """Carry out Fourier Domain Null Test PSF-corrected shape measurement routines.

    Example usage
    -------------

    Typical application to a single object:

        >>> galaxy = galsim.Gaussian(flux = 1.0, sigma = 1.0)
        >>> galaxy = galaxy.shear(g1=0.05, g2=0.0)  # shears the Gaussian by (0.05, 0) using the
        >>>                                         # |g| = (a - b)/(a + b) definition
        >>> galaxy_image = galaxy.drawImage(dx = 0.2)
        >>> result = fdnt.GLMoments(galaxy_image,)

    After running the above code, `result.observed_shape` ["shape" = distortion, the
    (a^2 - b^2)/(a^2 + b^2) definition of ellipticity] is
    `(0.XX, 0.XX)`  compared with the  expected `(0.09975, 0)` for a perfect measurement.

    @param gal_image       The Image of the galaxy being measured.
    @param weight          The optional weight image for the galaxy being measured.  Can be an int
                           or a float array.  Currently, GalSim does not account for the variation
                           in non-zero weights, i.e., a weight map is converted to an image with 0
                           and 1 for pixels that are not and are used.  Full use of spatial
                           variation in non-zero weights will be included in a future version of
                           the code.
    @param guess_x_wc      An initial guess for the x component of the object centroid (useful in
                           case it is not located at the center, which is the default
                           assumption).  The convention for centroids is such that the center of
                           the lower-left pixel is (0,0). [default: gal_image.trueCenter().x]
    @param guess_y_wc      An initial guess for the y component of the object centroid (useful in
                           case it is not located at the center, which is the default
                           assumption).  The convention for centroids is such that the center of
                           the lower-left pixel is (0,0). [default: gal_image.trueCenter().y]
    @param guess_sig_gal   Optional argument with an initial guess for the Gaussian sigma of the
                           galaxy (in pixels).
    @param order           The optional GL-fit order, when GL fits are used to fill masked-out
                           pixels.  [default: 0]
    @param bg              The optional background level. [default: 0.]
    @param sky             The optional sky level.  [default: 0.]

    @returns a ShapeData object containing the results of shape measurement.
    """

    if guess_a_b_pa == None and guess_g1g2 == None:
        raise RuntimeError("Either guess_a_b_pa or guess_g1g2 must be specified")
    elif guess_g1g2==None:
        (guess_a_wc, guess_b_wc, guess_pa_wc) = guess_a_b_pa
    else:  # assuming wc == pix, so guess_sig_gal_pix==guess_sig_gal_wc
        (guess_a_wc, guess_b_wc, guess_pa_wc) = calculate_a_b_pa(guess_sig_gal_pix, guess_g1g2)

    # prepare inputs to C++ routines:
    # _GLMoments(img::Image<float> gal_image, img::Image<float> weight_image,
    #            double x_pix, double y_pix, double a_wc, double b_wc, double pa_wc,
    #            double r_pix, double bg, int order, double sky)
    ## convert image formats
    gal_fdnt_image = _fdnt.FDNTImageF(gal_image.array, gal_image.xmin, gal_image.ymin)
    weight_fdnt_image = _convertMask(gal_image, weight=weight, badpix=badpix)
    ## convert int to float
    guess_x_wc = float(guess_x_wc)   # centroids often specified by integers (pixels)
    guess_y_wc = float(guess_y_wc)

    try:
        result = _fdnt._GLMoments(gal_fdnt_image, weight_fdnt_image,
                                  x_pix=guess_x_wc, y_pix=guess_y_wc,  # currently wc == pix ...
                                  a_wc=guess_a_wc, b_wc=guess_b_wc, pa_wc=guess_pa_wc,
                                  r_pix=guess_sig_gal_pix, bg=bg, order=order, sky=sky)

    except RuntimeError as err:
        raise RuntimeError

    return FDNTShapeData(result)


def calculate_a_b_pa(sigma, g1g2):
    """
    Estimate major axis a, minor axis b, and position angle pa (the SExtractor output format)
    The math goes as follows:
       (a^2-b^2)/(a^2+b^2) = |e|
       ab = sigma^2 (conservation of area after non-dilating shear)
       let   x == a/b > 1
       then  x^2 = (1+e)/(1-e)
             a^2 = sigma^2 * x
             b^2 = sigma^2 * x^-1
             pa = arctan(e_2, e_1) / 2 = arctan(g_2, g_1) / 2  (orientation is the same)
    """
    (g1, g2) = g1g2
    s = galsim.Shear(g1=g1, g2=g2)
    x_sqrt = (((1+s.getE()) / (1-s.getE()))) ** 0.25
    a = sigma * x_sqrt
    b = sigma / x_sqrt
    pa = (np.arctan2(g2, g1) / 2.0) / np.pi * 180.

    return a, b, pa
