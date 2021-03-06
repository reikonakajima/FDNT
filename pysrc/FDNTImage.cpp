
#include "boost/python.hpp"

#include "NumpyHelper.h"
#include "Image.h"

namespace bp = boost::python;
using img::ImageData;

namespace fdnt {


/*
 * Wrapper struct for ImageData<T>
 */
template <typename T>
struct PyFDNTImageData {

  static ImageData<T>* MakeFromArray(bp::object& array, int xmin, int ymin)
  {
    T* data = 0;                  // set in CheckNumpyArray()
    const int ndim = 2;           // array must be 2 dimensional
    bool isConst = true;          // won't work if false, but I don't understand why (ask Jim?)
    boost::shared_ptr<T> owner;   // ImageData is never the owner (someone else made the array)
    int stride = 0;               // set in CheckNumpyArray()
    CheckNumpyArray(array, ndim, isConst, data, owner, stride);

    int xsize = GetNumpyArrayDim(array.ptr(), 1);
    int ysize = GetNumpyArrayDim(array.ptr(), 0);
    Bounds<int> bounds(xmin, xmin+xsize-1, ymin, ymin+ysize-1);

    T *virtual_data_ptr = data - bounds.getXMin();    // yikes!  (inherited NR ugliness)
    T **true_row_ptrs = new T*[ysize];
    T **row_ptrs = true_row_ptrs - bounds.getYMin();  // yikes!
    for (int i=bounds.getYMin(); i<=bounds.getYMax(); i++) {
      row_ptrs[i] = virtual_data_ptr;
      virtual_data_ptr += xsize;
    }

    bool contiguous = true; // assuming GalSim generated image arrays will always be contiguous

    // Construct ImageData without owning the data (does not have array data delete previliges)
    return new ImageData<T>(bounds, row_ptrs, contiguous);
  }

  static bp::object GetArrayImpl(bp::object self, bool isConst)
  {
    // --- Try to get cached array ---
    if (PyObject_HasAttrString(self.ptr(), "_array") && self.attr("_array") != bp::object())
      return self.attr("_array");

    const ImageData<T>& imdata = bp::extract<const ImageData<T>&>(self);

    const T* data = imdata.location(imdata.getBounds().getXMin(), imdata.getBounds().getYMin());
    int n1 = imdata.getBounds().getYMax() - imdata.getBounds().getYMin() + 1;
    int n2 = imdata.getBounds().getXMax() - imdata.getBounds().getXMin() + 1;
    int stride = n2;  // this will NOT work for subimages!  TEMPORARY SOLUTION, checked by next line
    if (!imdata.contiguousData()) {
      throw ImageError("stride does not match X dimension");
    }
    bp::object numpy_array = MakeNumpyArray(data, n1, n2, stride, isConst);  // will always be owner.

    self.attr("_array") = numpy_array;
    return numpy_array;
  }

  static bp::object GetArray(bp::object imdata) { return GetArrayImpl(imdata, false); }



  static bp::object wrapFDNTImageData(const std::string& suffix) {

    typedef ImageData<T>* (*constructFromArray_func_type)(bp::object&, int, int);

    bp::class_< ImageData<T> >
      pyFDNTImageData(("FDNTImageData" + suffix).c_str(), "", bp::no_init);
    pyFDNTImageData
      .def(bp::init< Bounds<int>, T >(bp::args("bounds", "init_val")))
      .def(bp::init< Bounds<int> >(bp::args("bounds")))
      .def(
	   "__init__",
	   bp::make_constructor(
		constructFromArray_func_type(&MakeFromArray),
		bp::default_call_policies(), bp::args("array", "xmin", "ymin")
				)
	   )
      .add_property("array", &GetArray)  // getter only; set through ImageData constructor
      .def("getBounds", &ImageData<T>::getBounds)
      .add_property("bounds", &ImageData<T>::getBounds)
      ;

    // copy constructor from its own type
    //wrapFDNTImageDataTemplates<T>(pyFDNTImageData);

    return pyFDNTImageData;

  } // wrapFDNTImageData()

}; // struct PyFDNTImageData


/*
 * Wrapper struct for Image<T>
 */
template <typename T>
struct PyFDNTImage {

  static Image<T>* MakeFromArray(bp::object& array, int xmin, int ymin)
  {
    // Construct ImageData without owning the data (does not have data delete previliges)
    ImageData<T> *imdata_ptr = PyFDNTImageData<T>::MakeFromArray(array, xmin, ymin);

    ImageHeader *imhead_ptr = new ImageHeader();

    // Now construct the ImageData
    return new Image<T>(imdata_ptr, imhead_ptr);
  }

  static bp::object GetArrayImpl(bp::object self, bool isConst)
  {
    // --- Try to get cached array ---
    if (PyObject_HasAttrString(self.ptr(), "_array") && self.attr("_array") != bp::object())
      return self.attr("_array");

    const Image<T>& image = bp::extract<const Image<T>&>(self);

    const T* data = image.data()->location(image.getBounds().getXMin(), image.getBounds().getYMin());
    int n1 = image.getBounds().getYMax() - image.getBounds().getYMin() + 1;
    int n2 = image.getBounds().getXMax() - image.getBounds().getXMin() + 1;
    int stride = n2;  // this will NOT work for subimages!  TEMPORARY SOLUTION, checked by next line
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

    typedef Image<T>* (*constructFromArray_func_type)(bp::object&, int, int);
    typedef const T& (Image<T>::*at_func_type)(int, int) const;

    bp::object at = bp::make_function(
	at_func_type(&Image<T>::at),
	bp::return_value_policy<bp::copy_const_reference>(),
	bp::args("x", "y")
    );

    bp::class_< Image<T> >   // in Gary's code, it's an "Image<>" class, not "FDNTImage"
      pyFDNTImage(("FDNTImage" + suffix).c_str(), "", bp::no_init);
    pyFDNTImage
      .def(bp::init<int, int>(bp::args("ncol", "nrow")))
      .def(bp::init< Bounds<int>&, T >
	   ((bp::arg("bounds"), bp::arg("init_value")=T(0))))
      .add_property("array", &GetArray)  // getter only; set through ImageData constructor
      .def(
	   "__init__",
	   bp::make_constructor(
		constructFromArray_func_type(&MakeFromArray),
		bp::default_call_policies(), bp::args("array", "xmin", "ymin")
				)
	   )
      .def("getBounds", &Image<T>::getBounds)
      .add_property("bounds", &Image<T>::getBounds)
      .def("__call__", at) // always used checked accessors in Python
      ;

    // copy constructor from its own type
    wrapFDNTImageTemplates<T>(pyFDNTImage);

    return pyFDNTImage;

  } // wrapFDNTImage()

}; // struct PyFDNTImage



void pyExportFDNTImage() {

  bp::dict pyFDNTImageDataDict;

  pyFDNTImageDataDict[GetNumPyType<int16_t>()] = PyFDNTImageData<int16_t>::wrapFDNTImageData("S");
  pyFDNTImageDataDict[GetNumPyType<int32_t>()] = PyFDNTImageData<int32_t>::wrapFDNTImageData("I");
  pyFDNTImageDataDict[GetNumPyType<float>()] = PyFDNTImageData<float>::wrapFDNTImageData("F");
  pyFDNTImageDataDict[GetNumPyType<double>()] = PyFDNTImageData<double>::wrapFDNTImageData("D");

  bp::dict pyFDNTImageDict;

  pyFDNTImageDict[GetNumPyType<int16_t>()] = PyFDNTImage<int16_t>::wrapFDNTImage("S");
  pyFDNTImageDict[GetNumPyType<int32_t>()] = PyFDNTImage<int32_t>::wrapFDNTImage("I");
  pyFDNTImageDict[GetNumPyType<float>()] = PyFDNTImage<float>::wrapFDNTImage("F");
  pyFDNTImageDict[GetNumPyType<double>()] = PyFDNTImage<double>::wrapFDNTImage("D");

  bp::scope scope;  // a default constructed scope represents the module we're creating
  scope.attr("FDNTImageDataDict") = pyFDNTImageDataDict;
  scope.attr("FDNTImageDict") = pyFDNTImageDict;

}


}  // namespace fdnt
