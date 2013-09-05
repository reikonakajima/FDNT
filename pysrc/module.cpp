#include "boost/python.hpp"
#include "numpy/arrayobject.h"


// declare the export functions

namespace galsim {

  namespace fdnt {
    
    void pyExportRunFDNT();
    void pyExportFDNTImage();

  } 

}


// generate the galsim FDNT module

BOOST_PYTHON_MODULE(_fdnt) {

  galsim::fdnt::pyExportRunFDNT();
  galsim::fdnt::pyExportFDNTImage();

}
