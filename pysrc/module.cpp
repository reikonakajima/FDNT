#include "boost/python.hpp"

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#define PY_ARRAY_UNIQUE_SYMBOL FDNT_ARRAY_API
#include "numpy/arrayobject.h"


// declare the export functions

namespace fdnt {

  void pyExportBounds();
  void pyExportFDNTImage();
  void pyExportRunFDNT();


}


// generate the FDNT module

BOOST_PYTHON_MODULE(_fdnt) {

  import_array(); // for numpy
  fdnt::pyExportBounds();
  fdnt::pyExportRunFDNT();
  fdnt::pyExportFDNTImage();

}
