#include "ColorGrad.h"

ColorGradGaussian::ColorGradGaussian(string data_file_name) : psf() {

   // Read in PSF-info file
  ifstream psf_data(data_file_name.c_str());
  if (!psf_data) {
    cerr << "Error opening catalog file " << data_file_name << endl;
    exit(1);
  }
  string buffer;
  vector<double> lambda;
  vector<double> sigma;
  vector<double> norm;
  string l, s, n;
  while (stringstuff::getlineNoComment(psf_data, buffer)) {
    istringstream iss(buffer);
    if (!(iss >> l >> s >> n)) {
      cerr << "ColorGradGaussian: Error in line: " << buffer << endl;
      exit(1);
    }
    lambda.push_back(atof(l.c_str()));
    sigma.push_back(atof(s.c_str()));
    norm.push_back(atof(n.c_str()));
  }
  if (lambda.size() <= 0) {
    cerr << "ColorGradGaussian: no parameters read" << endl;
    exit(1);
  }

  // Convert to SBProfile sums
  for (int i = 0; i < lambda.size(); ++i) {
    SBGaussian sb_gauss(norm[i], sigma[i]);
    psf.add(sb_gauss);
  }
  cerr << "Size of ColorGradGaussian parameters .... " << lambda.size() << endl;
}
