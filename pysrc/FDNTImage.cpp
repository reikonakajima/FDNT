
#include "boost/python.hpp"

#include "NumpyHelper.h"
#include "Image.h"

namespace bp = boost::python;


namespace galsim {
namespace fdnt {

/*
 * Wrapper struct for Image Header
 */
struct PyFDNTImageHeader{



}; // struct PyFDNTImageHeader



/*
 * Wrapper struct for ImageData<T>
 */
template <typename T>
struct PyFDNTImageData {

  // XXX TODO: FILL ME XXX

}; // struct PyFDNTImageData


/*
 * Wrapper struct for Image<T>
 */
template <typename T>
struct PyFDNTImage {

  static bp::object GetArrayImpl(bp::object self, bool isConst)
  {
    // --- Try to get cached array ---
    if (PyObject_HasAttrString(self.ptr(), "_array") && self.attr("_array") != bp::object())
      return self.attr("_array");

    const Image<T>& image = bp::extract<const Image<T>&>(self);

    const T* data = image.data()->location(image.getBounds().getXMin(), image.getBounds().getYMin());
    int n1 = image.getBounds().getYMax() - image.getBounds().getYMin() + 1;
    int n2 = image.getBounds().getXMax() - image.getBounds().getXMin() + 1;
    int stride = n2;  // this will NOT work for subimages!  TEMOPRARY SOLUTION, checked by next line
    if (!image.data()->contiguousData()) {
      throw ImageError("stride does not match X dimension");
    }
    bp::object numpy_array = MakeNumpyArray(data, n1, n2, stride, isConst);  // will always be owner.

    self.attr("_array") = numpy_array;
    return numpy_array;
  }

  static bp::object GetArray(bp::object image) { return GetArrayImpl(image, false); }

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
      .add_property("array", &GetArray)
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
