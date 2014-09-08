// RunPSFEx.cpp
// Generate PSFEx psf models and measure shape; generate catalog

#include "StringStuff.h"
#include "FITSImage.h"
#include <fstream>
#include "PSFEx.h"
#include "Astrometry.h"
//#include "Image.h"
//#include "HeaderFromStream.h"

//using namespace laguerre;
using namespace sbp;
using namespace astrometry;

int
main(int argc, 
     char *argv[]) {
    
  // Input files
  string psfName = argv[1];
  string catName = argv[2];
  int psfOrder = -1;  // use PSFEx order
  
  try {

    // (1) Read the PSF
    PSFExModel *model;
    try {
      model = new PSFExModel(psfName.c_str());
    } catch (PSFExError &e) {
      cerr << "Error reading PSFEx model: " << e.what() << "; exiting." << endl; 
      exit(1);
    } catch (...) {
      cerr << "Error reading PSFEx model; exiting." << endl; 
      exit(1);
    }

    if (!model->selfTest(1, 3, psfOrder)) { // this is a very lenient self-test
      cerr << "PSFEx model did not pass self test; exiting." << endl;
      exit(1);
    }

    ifstream ccat(catName.c_str());
    if (!ccat) {
      cerr << "Error opening catalog file " << catName << endl;
      exit(1);
    }

    string buffer;
    while (stringstuff::getlineNoComment(ccat, buffer)) { // read x_pix, y_pix (position of galaxy)

      istringstream iss(buffer);
      double x_pix, y_pix;
      iss >> x_pix >> y_pix;
      if (!iss) {
	cerr << "Bad catalog input line: " << buffer;
	exit(1);
      }

      // generate model PSF at this location
      model->fieldPosition(x_pix, y_pix); // model now holds the psf model at this position
      cerr << "flux before normalization: " << model->sb()->getFlux() << endl;
      //model->setFlux(1.0);  // this is redundant
	    
      // generate PSF image in WC
      SBDistort psfWCS(*(model->sb()));
      psfWCS.setFlux(1.); 
	    
      double dx = 1.0;  // use ee50psf / 2.35 for fully sampled PSF
      cerr << "drawing psf postage stamp with dx=" << dx << endl;
      Image<> ipsf = psfWCS.draw(dx);
	    
      // draw onto a fits file for viewing

    }
    delete model;
	
  } catch (std::runtime_error &m) {
    cerr << m.what() << endl;
    quit(m,1);
  }
}
  
