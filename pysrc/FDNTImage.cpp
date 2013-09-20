
#include "boost/python.hpp"

#include "NumpyHelper.h"
#include "Image.h"

namespace bp = boost::python;


namespace galsim {
  namespace fdnt {


void pyExportFDNTImage() {

  // define Gary's old Image class with BoostPython
  // bp::class_<Image>("FDNTImage") ...;
  /*
  bp::class_<Image<float> >("FDNTImageF");
  bp::class_<Image<int> >("FDNTImageI");
  */
}


}  // namespace fdnt
}  // namespace galsim
