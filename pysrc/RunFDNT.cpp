
#include "boost/python.hpp"

#include "NumpyHelper.h"
#include "FDNT.h"
#include "Image.h"
#include "RunFDNT.h"

namespace bp = boost::python;
using namespace laguerre;
using namespace img;

namespace fdnt {
namespace {

struct PyFDNTShapeData {

  template <typename T>
  static void wrapTemplates() {
    typedef FDNTShapeData (*RunFDNT_func) (const Image<T>&, const Image<T>&, const Image<T>&,
					    double, double, double, double, double, double, double,
					    double, int, double);
    bp::def("_RunFDNT",
	    RunFDNT_func(&RunFDNT<T>),
	    (bp::arg("gal_image"), bp::arg("psf_image"), bp::arg("weight_image"),
	     bp::arg("x_pix"), bp::arg("y_pix"), bp::arg("a_wc"), bp::arg("b_wc"), bp::arg("pa_wc"),
	     bp::arg("r_pix"), bp::arg("ee50psf"), bp::arg("bg"), bp::arg("order"), bp::arg("sky")),
	    "Run the FDNT shape measurement method for a given galaxy and psf.");
  }  // static void wrapTemplates()


  static void wrap() {

    static char const * doc =
      "FDNTShapeData object represents information from the GLMoments and FDNT PSF-correction\n"
      "functions.  See C++ docs for more detail.\n"
      ;

    bp::class_<FDNTShapeData>("_FDNTShapeData", doc, bp::init<>())
      .def_readwrite("image_bounds", &FDNTShapeData::image_bounds)
      .def_readwrite("observed_flags", &FDNTShapeData::observed_flags)
      .def_readwrite("observed_e1", &FDNTShapeData::observed_e1)
      .def_readwrite("observed_e2", &FDNTShapeData::observed_e2)
      .def_readwrite("observed_sigma", &FDNTShapeData::observed_sigma)
      .def_readwrite("observed_b00", &FDNTShapeData::observed_b00)
      .def_readwrite("observed_b00_var", &FDNTShapeData::observed_b00_var)
      .def_readwrite("observed_b22", &FDNTShapeData::observed_b22)
      .def_readwrite("observed_centroid", &FDNTShapeData::observed_centroid)

      .def_readwrite("psf_flags", &FDNTShapeData::psf_flags)
      .def_readwrite("psf_e1", &FDNTShapeData::psf_e1)
      .def_readwrite("psf_e2", &FDNTShapeData::psf_e2)
      .def_readwrite("psf_sigma", &FDNTShapeData::psf_sigma)
      .def_readwrite("psf_order", &FDNTShapeData::psf_order)
      .def_readwrite("psf_b00", &FDNTShapeData::psf_b00)
      .def_readwrite("psf_b00_var", &FDNTShapeData::psf_b00_var)
      .def_readwrite("psf_b22", &FDNTShapeData::psf_b22)
      .def_readwrite("psf_chisq", &FDNTShapeData::psf_chisq)
      .def_readwrite("psf_DOF", &FDNTShapeData::psf_DOF)

      .def_readwrite("intrinsic_flags", &FDNTShapeData::intrinsic_flags)
      .def_readwrite("intrinsic_e1", &FDNTShapeData::intrinsic_e1)
      .def_readwrite("intrinsic_e2", &FDNTShapeData::intrinsic_e2)
      .def_readwrite("intrinsic_e1_var", &FDNTShapeData::intrinsic_e1_var)
      .def_readwrite("intrinsic_e2_var", &FDNTShapeData::intrinsic_e2_var)
      .def_readwrite("intrinsic_e1e2_covar", &FDNTShapeData::intrinsic_e1e2_covar)
      .def_readwrite("intrinsic_sigma", &FDNTShapeData::intrinsic_sigma)
      .def_readwrite("shrink_response", &FDNTShapeData::shrink_response)
      .def_readwrite("evaluation_count", &FDNTShapeData::evaluation_count)
      .def_readwrite("e_trial_count", &FDNTShapeData::e_trial_count)
      .def_readwrite("resolution_factor", &FDNTShapeData::resolution_factor)
      .def_readwrite("error_message", &FDNTShapeData::error_message)
      ;

    wrapTemplates<float>();
    //wrapTemplates<double>();
  }
};  // struct PyFDNTShapeData

} // namespace anonymous


void pyExportRunFDNT() {

  bp::def("greet", greet);
  PyFDNTShapeData::wrap();

}  // pyExportRunFDNT()


}  // namespace fdnt
