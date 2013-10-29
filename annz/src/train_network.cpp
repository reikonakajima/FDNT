#include "train_network.h"

//-----------------------------------------------------------------
// Cost function for individual training patterns.

/*
// Compute sum-of-squares error on given pattern.
double TrainNetwork::En(const std::vector<double>& inputs, const std::vector<double>& trues) {
    
  // Propagate pattern through the network.
  setInputs(inputs);
  updateActivations();
    
  // Sum squared error over outputs.
  double sumsq = 0.0;
  double dif;
  for (int i = 0; i < nOutputs; ++i) {
    dif = (nodes[outputLayer[i]]->activation - trues[i]);
    sumsq += dif * dif;
  }
    
  return sumsq;
}
*/

// Compute chi2 error on given pattern.
double TrainNetwork::En(const std::vector<double>& inputs, const std::vector<double>& trues, const std::vector<double>& sigma2s) {

  // Propagate pattern through the network.
  setInputs(inputs);
  updateActivations();

  // Sum squared error over outputs.
  double sumsq = 0.0;
  double dif;
  for (int i = 0; i < nOutputs; ++i) {
    dif = (nodes[outputLayer[i]]->activation - trues[i]);
    sumsq += dif * dif / sigma2s[i];
  }

  return sumsq;
}

/*
// Return derivative of cost function wrt weights, for the given training pattern.
std::vector<double> TrainNetwork::dEndw(const std::vector<double>& inputs, const std::vector<double>& trues) {
    
  // Must forward pass before updating deltas, in order to get activations.
  setInputs(inputs);
  updateActivations();
  updateDeltas(trues);
    
  std::vector<double> derivs;
    
  // Compute derivative wrt each weight in turn.
  // this is given by act_f * del_t
  double del,act;
  for (int i = 0; i < weights.size(); ++i) {
    del = nodes[weights[i]->to]->delta;
    act = nodes[weights[i]->from]->activation;
    derivs.push_back(del*act);
  }

  return derivs;
}
*/


// Return derivative of cost function wrt weights, for the given training pattern. This uses uncertainty of truth (chi^2 minimization)
std::vector<double> TrainNetwork::dEndw(const std::vector<double>& inputs, const std::vector<double>& trues, const std::vector<double>& sigma2s) {
    
  // Must forward pass before updating deltas, in order to get activations.
  setInputs(inputs);
  updateActivations();
  updateDeltas(trues,sigma2s);
    
  std::vector<double> derivs;
    
  // Compute derivative wrt each weight in turn.
  // this is given by act_f * del_t
  double del,act;
  for (int i = 0; i < weights.size(); ++i) {
    del = nodes[weights[i]->to]->delta;
    act = nodes[weights[i]->from]->activation;
    derivs.push_back(del*act);
  }

  return derivs;
}

