#
#
#
"""@file fdnt.py 
Routines for shape measurements with the FDNT method.

This file contains the python interface to C++ routines for shape measurement, given the galaxy,
PSF and mask images (reference: Bernstein 2010).

"""


from . import _fdnt
import fdntimage
import bounds

#from _fdnt import _FDNTParams as _FDNTParams

"""
class ShapeData(object):
    " ""A class to contain the outputs of using the FDNT shape and moments measurement routines.

    At the C++ level, we have a container for the outputs of the HSM shape measurement routines.
    The ShapeData class is the analogous object at the python level.  It contains the following
    information about moment measurement (from either EstimateShear() or FindAdaptiveMom()):

    - image_bounds: a BoundsI object describing the image.

    - moments_status: the status flag resulting from moments measurement; -1 indicates no attempt to
      measure, 0 indicates success.

    - observed_shape: a Shear object representing the observed shape based on adaptive
      moments.

    - moments_sigma: size sigma = (det M)^(1/4) from the adaptive moments, in units of pixels; -1 if
      not measured.

    - moments_amp: total image intensity for best-fit elliptical Gaussian from adaptive moments.
      Normally, this field is simply equal to the image flux (for objects that follow a Gaussian
      light distribution, otherwise it is something approximating the flux).  However, if the image
      was drawn using `drawImage(method='sb')` then moments_amp relates to the flux via
      flux = (moments_amp)*(pixel scale)^2.

    - moments_centroid: a PositionD object representing the centroid based on adaptive
      moments.  The convention for centroids is such that the center of the lower-left pixel is
      (0,0).

    - moments_rho4: the weighted radial fourth moment of the image.

    - moments_n_iter: number of iterations needed to get adaptive moments, or 0 if not measured.

    If EstimateShear() was used, then the following fields related to PSF-corrected shape will also
    be populated:

    - correction_status: the status flag resulting from PSF correction; -1 indicates no attempt to
      measure, 0 indicates success.

    - corrected_e1, corrected_e2, corrected_g1, corrected_g2: floats representing the estimated
      shear after removing the effects of the PSF.  Either e1/e2 or g1/g2 will differ from the
      default values of -10, with the choice of shape to use determined by the quantity meas_type (a
      string that equals either 'e' or 'g') or, equivalently, by the correction method (since the
      correction method determines what quantity is estimated, either the shear or the distortion).

    - corrected_shape_err: shape measurement uncertainty sigma_gamma per component.

    - correction_method: a string indicating the method of PSF correction (will be "None" if
      PSF-correction was not carried out).

    - resolution_factor: Resolution factor R_2;  0 indicates object is consistent with a PSF, 1
      indicates perfect resolution.

    - error_message: a string containing any error messages from the attempt to carry out
      PSF-correction.

    The ShapeData object can be initialized completely empty, or can be returned from the
    routines that measure object moments (FindAdaptiveMom()) and carry out PSF correction
    (EstimateShear()).
    " ""
    def __init__(self, *args):
        # arg checking: require either a CppShapeData, or nothing
        if len(args) > 1:
            raise TypeError("Too many arguments to initialize ShapeData!")
        elif len(args) == 1:
            if not isinstance(args[0], _galsim._CppShapeData):
                raise TypeError("Argument to initialize ShapeData must be a _CppShapeData!")
            self.image_bounds = args[0].image_bounds
            self.moments_status = args[0].moments_status
            self.observed_shape = galsim.Shear(args[0].observed_shape)
            self.moments_sigma = args[0].moments_sigma
            self.moments_amp = args[0].moments_amp
            self.moments_centroid = args[0].moments_centroid
            self.moments_rho4 = args[0].moments_rho4
            self.moments_n_iter = args[0].moments_n_iter
            self.correction_status = args[0].correction_status
            self.corrected_e1 = args[0].corrected_e1
            self.corrected_e2 = args[0].corrected_e2
            self.corrected_g1 = args[0].corrected_g1
            self.corrected_g2 = args[0].corrected_g2
            self.meas_type = args[0].meas_type
            self.corrected_shape_err = args[0].corrected_shape_err
            self.correction_method = args[0].correction_method
            self.resolution_factor = args[0].resolution_factor
            self.error_message = args[0].error_message
        else:
            self.image_bounds = _galsim.BoundsI()
            self.moments_status = -1
            self.observed_shape = galsim.Shear()
            self.moments_sigma = -1.0
            self.moments_amp = -1.0
            self.moments_centroid = _galsim.PositionD()
            self.moments_rho4 = -1.0
            self.moments_n_iter = 0
            self.correction_status = -1
            self.corrected_e1 = -10.
            self.corrected_e2 = -10.
            self.corrected_g1 = -10.
            self.corrected_g2 = -10.
            self.meas_type = "None"
            self.corrected_shape_err = -1.0
            self.correction_method = "None"
            self.resolution_factor = -1.0
            self.error_message = ""
"""

# A helper function for taking input weight and badpix Images, and returning a weight Image in the
# format that the C++ functions want
def _convertMask(image, weight = None, badpix = None):
    """Convert from input weight and badpix images to a single mask image needed by C++ functions.

    This is used by RunFDNT().
    """
    # if no weight image was supplied, make an int array (same size as gal image) filled with 1's
    if weight == None:
        b = _fdnt.BoundsI(image.bounds.xmin, image.bounds.xmax,  # convert types from
                          image.bounds.ymin, image.bounds.ymax)  # galsim.BoundsI to fdnt.BoundsI
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
        import numpy
        if weight.dtype == numpy.int32:
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


def RunFDNT(gal_image, PSF_image, guess_x_centroid, guess_y_centroid, 
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

    @param gal_image        The Image of the galaxy being measured.
    @param PSF_image        The Image for the PSF.
    @param weight           The optional weight image for the galaxy being measured.  Can be an int
                            or a float array.  Currently, GalSim does not account for the variation
                            in non-zero weights, i.e., a weight map is converted to an image with 0
                            and 1 for pixels that are not and are used.  Full use of spatial
                            variation in non-zero weights will be included in a future version of
                            the code.
    @param guess_x_centroid  An initial guess for the x component of the object centroid (useful in
                            case it is not located at the center, which is the default
                            assumption).  The convention for centroids is such that the center of
                            the lower-left pixel is (0,0). [default: gal_image.trueCenter().x]
    @param guess_y_centroid  An initial guess for the y component of the object centroid (useful in
                            case it is not located at the center, which is the default
                            assumption).  The convention for centroids is such that the center of
                            the lower-left pixel is (0,0). [default: gal_image.trueCenter().y]
    @param guess_sig_gal    Optional argument with an initial guess for the Gaussian sigma of the
                            galaxy (in pixels).
    @param guess_sig_PSF    Optional argument with an initial guess for the Gaussian sigma of the
                            PSF (in pixels).
    @param order            The optional GL-fit order, when GL fits are used to fill masked-out
                            pixels.  [default: 0]
    @param bg               The optional background level. [default: 0.]
    @param sky              The optional sky level.  [default: 0.]

    @returns a ShapeData object containing the results of shape measurement.
    """
    # prepare inputs to C++ routines: ImageView for galaxy, PSF, and weight map
    gal_fdnt_image = _fdnt.FDNTImageF(gal_image.array, gal_image.xmin, gal_image.ymin)
    PSF_fdnt_image = _fdnt.FDNTImageF(PSF_image.array, PSF_image.xmin, PSF_image.ymin)
    weight_fdnt_image = _convertMask(gal_image, weight=weight, badpix=badpix)

    try:
        result = _fdnt._RunFDNT(gal_fdnt_image, PSF_fdnt_image, weight_fdnt_image,
                                x_pix=guess_x_centroid, y_pix=guess_y_centroid,
                                r_pix=guess_sig_gal_pix, ee50psf=guess_sig_PSF_pix,
                                a_wc=guess_a_wc, b_wc=guess_b_wc, pa_wc=guess_pa_wc,
                                bg=bg, order=order, sky=sky)
                               
    except RuntimeError as err:
        raise RuntimeError

    return result  # should be 0 if successful



# make FindAdaptiveMom a method of Image class
##galsim.Image.FindAdaptiveMom = FindAdaptiveMom
