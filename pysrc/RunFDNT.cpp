
#include "boost/python.hpp"

#include "NumpyHelper.h"
#include "FDNT.h"
#include "Image.h"
#include "RunFDNT.h"

namespace bp = boost::python;


namespace galsim {
  namespace fdnt {


void pyExportRunFDNT() {

  // define the run code of FDNT
  //
  bp::def("greet", greet);

}


}  // namespace fdnt
}  // namespace galsim
