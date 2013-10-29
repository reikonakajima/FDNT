#include "network.h"
#include <cstdlib>

//----------------------------------------------------------------
// Network wrapper class, for use when training a new network.
//----------------------------------------------------------------
class TrainNetwork : public Network {
  
public:

  // Constructor.
  TrainNetwork(const char arch_file[]) : Network(arch_file) {
    // Assign random initial values to the weights.
    for (int w = 0; w < weights.size(); ++w) {
      weights[w]->wt = 2.0 * rand()/RAND_MAX - 1.0; 
    }
  }

  // Return En for given training pattern. (sum of squared errors)
  double En(const std::vector<double>& inputs, const std::vector<double> & trues, const std::vector<double> & sigma2s);
  
  // Return derivative of cost function wrt weights, for the given training pattern.
  std::vector<double> dEndw(const std::vector<double> & inputs, const std::vector<double> & trues, const std::vector<double> & sigma2s);

  // Change weight values.
  void setWeights(std::vector<double> wts) {
    for (int i = 0; i < weights.size(); ++i)
      weights[i]->wt = wts[i];
  }
  
  // Retrieve current weights vector.
  std::vector<double> getWeights() {
    std::vector<double> wts;
    for (int i = 0; i < weights.size(); ++i) {
      wts.push_back(weights[i]->wt);
    }
    return wts;
  }

};
