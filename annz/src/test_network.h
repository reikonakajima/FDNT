#include "network.h"

//-----------------------------------------------------------------------
// Network wrapper class, for use when applying network to test data.
//-----------------------------------------------------------------------
class TestNetwork : public Network {

public:
  
  // Constructor.
  TestNetwork(const char wts_file[]) : Network(wts_file) {

    // Attempt to restore saved weights and normalisation table.
    // (The call to the super constructor will have already read the architecture)
    restore_wts(wts_file);

    // We can now initialise the Jacobian std::vector of each node
    // to have dimension nInputs.
    for (int n = 0; n < nodes.size(); ++n) {
      nodes[n]->jacob.assign(nInputs, 0.0);
    }
  }
  
  // Restore saved weight values from the supplied file.
  bool restore_wts(const char wts_file[]);
  
  // Normalise the given input vector; pass through the network; 
  // renormalise the outputs; return results in the supplied vector.
  void norm_pass(std::vector<double> inputs, std::vector<double> & outputs, 
		 std::vector< std::vector<double> > & jacobian);
  
  // Normalisations.
  std::vector<double> meanIn, sigmaIn;
  std::vector<double> meanOut, sigmaOut;
  
};
