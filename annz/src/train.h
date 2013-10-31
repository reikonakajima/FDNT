#ifndef TRAIN_H
#define TRAIN_H

#include "train_network.h"
#include "util.h"

#include <vector>
#include <iostream>

// This class controls the training process.
class Training {

  // The network.
  TrainNetwork myNetwork;
  int nInputs, nOutputs; // Architecture.
  
  // Training pattern struct.
  struct Pattern {
    std::vector<double> inputs;
    std::vector<double> trues;
    std::vector<double> sigma2s; // sigma^2 uncertainty of truth
  };

  // The training and validation sets.
  std::vector<Pattern> trainSet, validSet;
  
  // Read in dataset from given filename.
  static bool read_dataset(char fileName[], int nIn, int nOut, std::vector<Pattern> & data, double sigma2tolerance);
  
  // Normalisation parameters.
  std::vector<double> meanIn, sigmaIn;
  std::vector<double> meanOut, sigmaOut;
  
  // Determine normalisation parameters and apply to training and validation sets.
  void normalise();

  // Cost function 'hyperparameters': the combination alpha/beta gives the weight decay parameter,
  // and it is this value which we report to the user.
  double alpha, beta;

  // Return cost function evaluated on training set, for the given set of weights.
  // This is the function to be minimised during training.
  double E(const std::vector<double> & wts);
  
  // Compute the derivatives of the cost function wrt the weights.
  // Returned in dEdw.
  void dEdw(const std::vector<double> & wts, std::vector<double> & dEdw);
  
  // Variable metric minimisation routine.
  void vmmin(int iter);

  // Re-evaluate hyperparameters (i.e. weight decay factor) at current training state.
  void recompute_hyperparams(const std::vector<double> & wts, std::vector< std::vector<double> > & G);

  // Used to store the best weight vector (that which gives minimum RMS on the validation set).
  std::vector<double> best_wts;
  double val_best; // RMS error for these best weights.

public:
  
  // Constructor.
  Training(char archFile[], char trainFile[], char validFile[], double sigma2tolerance);
  
  // Conduct nIter iterations of training.
  int train(int nIter) { vmmin(nIter); }
  
  // Compute sum-of-squares error on given set (or RMS error).
  // N.B. NOT THE SAME AS THE COST FUNCTION: NO WEIGHT DECAY.
  double sum_of_squares(const std::vector<Pattern> & data);
  double rms_error(const std::vector<Pattern> & data);
  
  // Save current weight values to file.
  void saveNetState(char file_name[]);

};

#endif // TRAIN_H

