//-----------------------------------------------------------------------------
// Implementation of the Artificial Neural Network.
//------------------------------------------------------------------------------

#ifndef NETWORK_H
#define NETWORK_H

#include <vector>
#include <cmath>
#include <iostream>

//----------------------------------------------------------------
// This class implements the nuts and bolts of the ANN.
// Although it may be instantiated, derived classes which provide 
// an extra layer of abstraction (optimised respectively for 
// training and testing a network) are provided in train_network
// test_network
//----------------------------------------------------------------
class Network {

private:
  
  // Node class.
  class Node {
  public:
    double activation;  // Activation of this node (Bishop's 'z').
    double delta;       // Backpropagation error.
    std::vector<double> jacob; // Jacobian for this node (derivatives wrt network inputs).
    Node() : activation(0), delta(0) {}           // Default constructor.
    Node(double a) : activation(a), delta(0) {}   // Constructor with initial activation.
  }; 

  // Connection type: wraps a weight value with the IDs of the nodes it connects.
  struct Connection {
    Connection(int t, int f)
      : from(f), to(t), wt(0.0) {}
    int from, to;
    double wt;
  };

  // Logistic activation function.
  // Input is activity (Bishop's 'a').
  static double logistic(const double a) {
    return 1.0/(1.0 + std::exp(-a));
  }

protected:

  // Network Layers - a "Layer" object contains a list of Node IDs.
  typedef std::vector<int> Layer;
  Layer inputLayer;
  Layer outputLayer;
  std::vector<Layer> hiddenLayer;
    
  // This network's nodes, indexed by their ID.
  std::vector<Node*> nodes;
  
  const static int BIAS_ID = 0; // The bias node always appears as the first entry in the nodes vector.

  // The weights vector.
  std::vector<Connection*> weights;

  // Lookup table for weights vector, indexed by (to, from).
  std::vector< std::vector<int> > wt_lu;

  // Methods used by Constructor to restore a Network from specification in file.
  bool restore_arch(const char file_name[]);
  void setup(int nIn, int nOut, int nLay, std::vector<int> nHid);
  
  // Note that a limited set of architectures are permitted: adjacent
  // layers are fully connected, and no connections are allowed between 
  // non-adjacent layers. Bishop makes it clear that such networks retain 
  // universality. Under these restrictions, the following are sufficient 
  // to completely specify the architecture.
  int nInputs; // Number of inputs.
  int nOutputs; // Number of outputs.
  int nLayers; // Number of hidden layers.
  std::vector<int> nHidden; // Number of nodes in each hidden layer.
  
  // Set inputs to values given in the vector.
  void setInputs(const std::vector<double>& inputs) {
    for (int i = 0; i < nInputs; ++i)
      nodes[inputLayer[i]]->activation = inputs[i];
  }
  
  // Return current activation of output nodes.
  std::vector<double> getOutputs() {
    std::vector<double> outs;
    for (int n = 0; n < outputLayer.size(); ++n)
      outs.push_back(nodes[outputLayer[n]]->activation);
    return outs;
  }
  
  // Update activations of all nodes, e.g. to reflect new inputs/weights.
  void updateActivations();
  
  // Update Jacobian at each node.
  void updateJacobian();
  
  // Update backpropagation errors.
  void updateDeltas(const std::vector<double> & trues, const std::vector<double> & sigma2s);

  // Constructor: restores network architecture and weight values (if present)
  // from the given file.  
  Network(const char file_name[]);

public:  
  
  // Set the network inputs to the given values, and propagate through the network,
  // returning the results in 'outputs'.
  void pass(const std::vector<double> & inputs, std::vector<double> & outputs);
  
  // Accessors for architecture.
  int getNInputs() const { return nInputs; }
  int getNOutputs() const { return nOutputs; }
  int getNWeights() const { return weights.size(); }
  
  // Dump weights to given output stream.
  void dump_weights(std::ostream& out);
  
  // Dump network architecture to given output stream.
  void dump_arch(std::ostream& out);

};

#endif // NETWORK_H
