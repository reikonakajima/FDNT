#include <iostream>
#include <string>
#include <sstream>
#include <fstream>

using namespace std;

int main()
{
  cout << "==============================" << endl;
  cout << " ANNz: Network architecture" << endl;
  cout << "==============================" << endl;

  //----------------------------------------------------
  // Get architecture specification from user
  //----------------------------------------------------
 
  int n_filt = 0;
  do {
    cout << "Number of filters (minimum 1): ";
    cin >> n_filt;
    cout << endl;
  } while (n_filt < 1);
  
  int n_hidlayer = 0;
  do {
    cout << "Number of hidden layers (minimum 1): ";
    cin >> n_hidlayer;
    cout << endl;
  } while (n_hidlayer < 1);

  // Get number of nodes in each hidden layer in turn
  int n_hidden[n_hidlayer];
  for (int ii = 0; ii < n_hidlayer; ++ii) {
    do {
      cout << "Nodes in hidden layer " << ii + 1 << ": ";
      cin >> n_hidden[ii];
      cout << endl;
    } while (n_hidden[ii] < 1);
  }

  // Number of outputs
  int n_outputs = 0;
  do {
    cout << "Number of outputs (minimum 1): ";
    cin >> n_outputs;
    cout << endl;
  } while (n_outputs < 1);

  //-------------------------------------
  // Construct architecture filename.
  //-------------------------------------
  stringstream ss(stringstream::out | stringstream::in);

  ss << "arch." << n_filt << ".";
  for (int ii = 0; ii < n_hidlayer; ++ii) {
    ss << n_hidden[ii] << ".";
  }
  ss << n_outputs << ".net";

  // Copy from stringstream to char*.
  const char * arch = ss.str().c_str();

  //---------------------------------------------
  // Write architecture to file.
  //---------------------------------------------
  
  ofstream arch_file(arch);
  
  arch_file << "# Generated by annz_net" << endl;
  arch_file << "# " << arch << endl;
  arch_file << "N_INPUTS " << n_filt << endl;
  arch_file << "N_OUTPUTS " << n_outputs << endl;
  arch_file << "N_LAYERS " << n_hidlayer << " ";
  for (int i = 0; i < n_hidlayer; ++i) {
    arch_file << n_hidden[i] << " ";
  }
  arch_file << endl;
  arch_file.close();

  cout << "======================================" << endl;
  cout << "Created \"" << arch << "\"" << endl;
  cout << "======================================" << endl;

  
}
