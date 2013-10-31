#include "train.h"

#include <iostream>

using namespace std;

int main(int argc, char * argv[]) {
  
  cerr.imbue(locale(""));


  cerr << endl << "===================================" << endl;
  cerr << " ANNz: Network training " << endl;
  cerr << "===================================" << endl << endl;
  
  // Check for minimum number of arguments
  if (argc < 7) {
    cerr << "ERROR: Insufficient arguments" << endl;
    cerr << "Usage: annz_train <arch_file> <train_file> <valid_file> <out_file> <rand_seed> <tolerance rms>" << endl;
    return -1;
  }

  // Seed random number generator (used to initialise the weights).
  srand(atol(argv[5]));

  // Create Training object.
  Training myTraining(argv[1], argv[2], argv[3], pow(atof(argv[6]),2));
  
  // Loop until no further iterations are requested.
  int niter;
  do {
    cout << "Maximum iterations: ";
    cin >> niter;
    cout << endl;
    if (niter == 0) break;
    myTraining.train(niter);
  } while (true);
  
  // Save best weight values to file.
  myTraining.saveNetState(argv[4]);

}
