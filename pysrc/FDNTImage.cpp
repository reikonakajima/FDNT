
#include "boost/python.hpp"

#include "NumpyHelper.h"
#include "Image.h"

namespace bp = boost::python;


namespace galsim {
  namespace fdnt {


void pyExportFDNTImage() {

  // define Gary's old Image class with BoostPython
  // bp::class_<Image>("FDNTImage") ...;
}


}  // namespace fdnt
}  // namespace galsim
