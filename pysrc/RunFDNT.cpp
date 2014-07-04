
#include "boost/python.hpp"

#include "NumpyHelper.h"
#include "FDNT.h"
#include "Image.h"
#include "RunFDNT.h"

namespace bp = boost::python;
using namespace laguerre;
using namespace img;

namespace fdnt {


struct PyRunFDNT {

  template <typename T>
  static void wrapTemplates() {
    typedef int (*RunFDNT_func) (const Image<T>&, const Image<T>&, const Image<T>&,
				 double, double, double, double, double, double, double,
				 double, int, double);
    bp::def("_RunFDNT", RunFDNT_func(&RunFDNT<T>),
	    (bp::arg("gal_image"), bp::arg("psf_image"), bp::arg("weight_image"),
	     bp::arg("x_pix"), bp::arg("y_pix"), bp::arg("a_wc"), bp::arg("b_wc"), bp::arg("pa_wc"),
	     bp::arg("r_pix"), bp::arg("ee50psf"), bp::arg("bg"), bp::arg("order"), bp::arg("sky")),
	    "Run the FDNT shape measurement method for a given galaxy and psf."
	    );
  }  // static void wrapTemplates()

};  // struct PyRunFDNT


void pyExportRunFDNT() {

  bp::def("greet", greet);
  PyRunFDNT::wrapTemplates<float>();
  /*
  PyRunFDNT::wrapTemplates<double>();
  PyRunFDNT::wrapTemplates<int>();
  PyRunFDNT::wrapTemplates<short>();
  */

}  // pyExportRunFDNT()


}  // namespace fdnt
