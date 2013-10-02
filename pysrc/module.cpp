#include "boost/python.hpp"

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#define PY_ARRAY_UNIQUE_SYMBOL FDNT_ARRAY_API
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
  import_array(); // for numpy
  galsim::fdnt::pyExportRunFDNT();
  galsim::fdnt::pyExportFDNTImage();

}
