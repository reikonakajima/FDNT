
#include "boost/python.hpp"

#include "NumpyHelper.h"
#include "Image.h"

namespace bp = boost::python;


namespace galsim {
namespace fdnt {

/*
 * Wrapper struct for Image<T>
 */
template <typename T>
struct PyFDNTImage {

  static Image<T>* MakeFromImage(const Image<T>& rhs)
  { return new Image<T>(rhs); }

  template <typename U, typename W>
  static void wrapFDNTImageTemplates(W& wrapper) {
    typedef Image<T>* (*constructFrom_func_type)(const Image<U>&);
    //typedef void (<Image<T>::* copyFrom_func_type)(const Image<U>&);
    wrapper
      .def(
	   "__init__",
	   bp::make_constructor(
		constructFrom_func_type(&MakeFromImage),
		bp::default_call_policies(), bp::args("other")
				)
	   )
      //.def("copyFrom", copyFrom_func_type(&Image<T>::copyFrom));
      ;
  } // wrapFDNTImageTemplates()

  static bp::object wrapFDNTImage(const std::string& suffix) {

    bp::class_< Image<T> >   // in Gary's code, it's an "Image<>" class
      pyFDNTImage(("FDNTImage" + suffix).c_str(), "", bp::no_init);
    pyFDNTImage
      .def(bp::init<int, int>(bp::args("ncol", "nrow")))
      ;

    // These lines allows for copy constructors between different types.  None allowed for now.
    //wrapFDNTImageTemplates<int16_t>(pyFDNTImage);
    //wrapFDNTImageTemplates<int32_t>(pyFDNTImage);
    //wrapFDNTImageTemplates<float>(pyFDNTImage);
    //wrapFDNTImageTemplates<double>(pyFDNTImage);

    // copy constructor from its own type
    wrapFDNTImageTemplates<T>(pyFDNTImage);

    return pyFDNTImage;

  } // wrapFDNTImage()

}; // struct PyFDNTImage



void pyExportFDNTImage() {

  /*
  PyFDNTImageHeader::wrapFDNTImageHeader;  // returns "FDNTImageHeader" class

  PyFDNTImageData<int16_t>::wrapFDNTImageData("S");  // returns "FDNTImageDataS" class
  PyFDNTImageData<int32_t>::wrapFDNTImageData("I");  // returns "FDNTImageDataI" class ... etc
  PyFDNTImageData<float>::wrapFDNTImageData("F");
  PyFDNTImageData<double>::wrapFDNTImageData("D");
  */

  PyFDNTImage<int16_t>::wrapFDNTImage("S");  // returns "FDNTImageS" class
  PyFDNTImage<int32_t>::wrapFDNTImage("I");  // returns "FDNTImageI" class ... etc
  PyFDNTImage<float>::wrapFDNTImage("F");
  PyFDNTImage<double>::wrapFDNTImage("D");

}


}  // namespace fdnt
}  // namespace galsim
