#include "train.h"

#include <fstream>
#include <iomanip>

using namespace std;

//============================================================================
// Training process.
//============================================================================

//--------------------------------------------------------------------------------
// The variable metric minimisation routine.
// Adapted from B.D. Ripley's C code.
void Training::vmmin(const int iter) {
  
  std::vector<double> g, gl, w, wl; // Weights and gradients.
  w = myNetwork.getWeights(); // Get initial weights values.
  int nwts = w.size();

  // Inverse Hessian approximation (nwts x nwts matrix).
  std::vector< std::vector<double> > G(nwts, std::vector<double>(nwts, 0.0));
  
  // Current value and previous minimum of the cost function.
  double f = E(w);
  double Fmin = f;

  std::cerr << "Initial cost function: " << std::fixed << std::setprecision(3) << f;

  dEdw(w, g); // Get initial gradients.

  int ilast = 0, hlast = 0, gradcount = 0;
  
  // Projection of the quasi-Newton direction on gradient vector
  double gradproj;
  std::vector<double> t(nwts, 0.0); // Quasi-Newton direction vector.

  double s; // Utility double for making sums.

  // Validation set error.
  double val_rms = rms_error(validSet);
  
  // Display validation set error (multiplying by the normalisation makes it easy to understand in terms
  // of rms photo-z error, but doesn't really make sense if we have multiple outputs).
  std::cerr << "  reduced chi^2 (valid):" << std::fixed << std::setprecision(5) << std::setw(10) << val_rms << std::endl;
  
  // Flags...
  bool accpoint, enough, reset = true, converged = false;

  // Parameters.
  const static double acctol = 1e-4, stepredn = 0.2, reltest = 10.0, abstol = 1e-2, reltol = 1e-8;
  
  //-------------------------------------------------------------------------------------
  // Main training loop.
  do {

    //-----------------------------------------------------
    // Check if reset is required.
    if (reset) {
      reset = false;
      ilast = gradcount;
      // Reset inverse Hessian to unit matrix -- 
      // this causes the next step to be taken in the simple steepest descent direction.
      for (int i = 0; i < nwts; i++) {
	// Only ever use lower triangle.
	for (int j = 0; j < i; j++)
	  G[i][j] = 0.0;
	G[i][i] = 1.0;
      }
      if (gradcount > 1) 
	std::cerr << "restarting" << std::endl;
    }

    // Copy current weights and derivatives.
    gl = g; wl = w;

    //----------------------------------------------------------------------
    // Compute the quasi-Newton direction -Gg and its projection onto gradient
    gradproj = 0.0;
    for (int i = 0; i < nwts; ++i) {
      s = 0.0;
      // Two-stage loop required because only lower triangle of G is updated.
      for (int j = 0; j <= i; j++)
	s -= G[i][j] * g[j];
      for (int j = i+1; j < nwts; j++)
	s -= G[j][i] * g[j];
      t[i] = s;
      gradproj += s * g[i];
    }

    //------------------------------------------------------------------
    // Check for quasi-Newton direction going uphill.
    if (gradproj >= 0.0) { 
      std::cerr << "Search uphill -- ";

      if (ilast == gradcount) {
	// We have just tried resetting.
	// This means that even a simple gradient descent cannot find a 
	// downhill direction -> in a minimum.
	converged = true;
	break;
      } else {
	reset = true;  // Try a reset
	continue;
      }      
    }

    //------------------------------------------------------------
    // In downhill direction - determine good step length.
    
    // Full Newton step: may be too large if quasi-Newton direction is wayward
    double steplength = 1.0; 
    
    accpoint = false; // Set to true once satisfactory reduction in cost function is obtained

    // Loop, reducing steplength until either the cost function is sufficiently reduced
    // or steplength becomes so small the weights do not change significantly.
    do {
	
      int count = 0;
	  
      // Take step in weight space
      for (int i = 0; i < nwts; i++) {
	w[i] = wl[i] + steplength * t[i]; // Update weights using current steplength
	if (reltest + w[i] == reltest + wl[i]) // Check for no significant change
	  count++; // Count how many weights have not changed
      }
	  
      if (count < nwts) { 
	// We've moved in weight space - have we reduced the cost function?

	f = E(w); // Re-evaluate cost function
	    
	// (gradproj * steplength) is the reduction in f we would expect if the error surface were linear
	// Minimum required reduction is specified by acctol
	accpoint = (f <= Fmin + gradproj * steplength * acctol);

	if (!accpoint) {
	  // f not sufficiently reduced -- reduce step length
	  steplength *= stepredn;
	}
      } else {
	// No significant movement in weight space.
	reset = true;
      }
	
      // Stop if the steplength is so small we're no longer moving significantly,
      // or an acceptable reduction in cost function is achieved  
    } while (!(reset || accpoint));
    
    // Force reset if value is small (in absolute terms) or if relative change is low
    enough = (f > abstol) && (f < (1.0-reltol)*(Fmin));
    if (!enough) reset = true;

    //---------------------------------------------------------------
    // Check whether reset has been requested (due to no progress).
    if (reset) {
      std::cerr << "Apparent convergence -- ";
      if (ilast == gradcount) {
	converged = true;
	break; // We just tried reset, with no success.
      }
    }
 
    //-----------------------------------------------------
    // We've made progress.
    
    Fmin = f; // This is now our minimum value
    gradcount++; // Update iteration counter.

    //-------------------------------------------------------
    // Report progress.
    std::cerr << "Iteration:" << std::setw(6) << gradcount;
    std::cerr << "     Train cost function:" << std::right << std::fixed << std::setw(12) << std::setprecision(3) << Fmin;
    
    val_rms = rms_error(validSet); // Get validation set error.
    std::cerr << "     reduced chi^2 (valid):" << std::fixed << std::setprecision(5) << std::setw(10) << val_rms;
    if (val_rms < val_best) {
      val_best = val_rms;
      best_wts = w;
      std::cerr << "  (best)";
    }
    std::cerr << std::endl;
      
    
    //-----------------------------------------------------------
    // Refine inverse Hessian approximation using BFGS update formula.
    
    dEdw(w, g); // Update gradient
    
    double D1 = 0.0;

    for (int i = 0; i < nwts; i++) {
      t[i] = steplength * t[i]; // My 't' == Bishop's 'p' -- i.e. the change in the weights
      gl[i] = g[i] - gl[i];
      D1 += t[i] * gl[i];
    }

    if (D1 > 0) {
	  
      double D2 = 0.0;
      for (int i = 0; i < nwts; i++) {
	s = 0.0;
	for (int j = 0; j <= i; j++)
	  s += G[i][j] * gl[j];
	for (int j = i+1; j < nwts; j++)
	  s += G[j][i] * gl[j];
	wl[i] = s;
	D2 += s * gl[i];
      }
	  
      D2 = 1.0 + D2 / D1;
      
      for (int i = 0; i < nwts; i++) {
	for (int j = 0; j <= i; j++)
	  G[i][j] += (D2 * t[i] * t[j] - wl[i] * t[j] - t[i] * wl[j]) / D1;
      }
	  
    } else { // D1 < 0
      std::cerr << "Not well-conditioned -- ";
      reset = true;
    }
    
    if (gradcount - ilast > 2*nwts) { reset = true; std::cerr << "Periodic reset -- ";}
    if (gradcount - hlast > 29) {
      hlast = gradcount;
      // Update alpha and beta.
      recompute_hyperparams(w, G);
      // The new hyperparameters will mean the cost function has changed.
      // re-evaluate, to avoid confusing the minimiser (if the cost function suddenly increases
      // it may claim convergence).
      Fmin = E(w);
    }

  } while (gradcount < iter);
  
  if (converged)
    std::cerr << "converged." << std::endl;
}

//----------------------------------------------------------
// RMS error on data set. Actually, this is the reduced chi^2
double Training::rms_error(const std::vector<Pattern> & data) {
  return sum_of_squares(data)/data.size();
}


//-------------------------------------------------------------
// Returns sum of squares error on the given data set.
double Training::sum_of_squares(const std::vector<Pattern> & data) {
    
  double total_sumsq = 0.0;
  
  // Compute squared error for each training pattern.
  for (int i = 0; i < data.size(); ++i) { 
    total_sumsq += myNetwork.En(data[i].inputs, data[i].trues, data[i].sigma2s); // En is hacked to be chi^2
  }

  return total_sumsq; 
}



//-------------------------------------------------------------------------------
// Return cost function evaluated on training set, for the given set of weights.
// This is the function to be minimised.
double Training::E(const std::vector<double> & wts) {

  myNetwork.setWeights(wts);

  // Calculate sum of squares on training set.
  double sumsq = sum_of_squares(trainSet);

  // Calculate weight decay.
  double wtDecay = 0.0;
  for (int i = 0; i < wts.size(); ++i) {
    wtDecay += wts[i] * wts[i];
  }

  //cout << "sumsq=" << beta * 0.5 * sumsq << " wtDecay=" << alpha * 0.5 * wtDecay << endl;
  
  // Return total cost function.
  return beta * 0.5 * sumsq + alpha * 0.5 * wtDecay;
}


//-----------------------------------------------------------------
// Computes the derivatives of the cost function wrt the weights.
// Returned in dEdw.
void Training::dEdw(const std::vector<double> & wts, std::vector<double> & dEdw) {

  myNetwork.setWeights(wts);

  std::vector<double> derivs(wts.size(), 0.0);
  std::vector<double> temp;
  
  for (int i = 0; i < trainSet.size(); ++i) {
    temp = myNetwork.dEndw(trainSet[i].inputs, trainSet[i].trues, trainSet[i].sigma2s); // dEndw is hacked to differentiate chi^2
    for (int j = 0; j < derivs.size(); ++j) {
      derivs[j] += temp[j];
    }
  }

  dEdw.assign(wts.size(), 0.0);
  for (int i = 0; i < wts.size(); ++i) {
    dEdw[i] = 0.5 * beta * derivs[i] + alpha * wts[i];
  }

}

//-----------------------------------------------------------------
// Update hyperparameters.
void Training::recompute_hyperparams(const std::vector<double> & wts, std::vector< std::vector<double> > & G) {
  
  cerr << "Updating weight decay parameter... ";

//   double gamma = 0;

//   // Evaluate trace of inverse Hessian.
//   for (int i = 0; i < G.size(); ++i) {
//     gamma += G[i][i];
//   }

//   gamma = G.size() - alpha * gamma;

//   if (gamma < 0.0) gamma = 1;

//   cerr << gamma << " of " << G.size() << " weights well-determined; ";
  
//   beta = (trainSet.size() * nOutputs - gamma)/sum_of_squares(trainSet);
  
  // Calculate weight decay.
  double wtDecay = 0.0;
  
  for (int i = 0; i < wts.size(); ++i)
    wtDecay += wts[i] * wts[i];
  
  //alpha = gamma/wtDecay;

  // "Quick and dirty"...
  alpha = wts.size()/wtDecay;

  beta = trainSet.size()/sum_of_squares(trainSet);

  cerr << "beta = " << alpha/beta << endl;

}

//==========================================================================
// Training setup.
//==========================================================================

// Constructor.
Training::Training(char archFile[], char trainFile[], char validFile[], double sigma2tolerance) :
      myNetwork(TrainNetwork(archFile)), alpha(1), beta(1) {

  // Set local record of numbers of inputs and outputs.
  nInputs = myNetwork.getNInputs();
  nOutputs = myNetwork.getNOutputs();
  
  // Read training set.
  read_dataset(trainFile, nInputs, nOutputs, trainSet, sigma2tolerance);
  std::cerr << "Training set: " << trainFile << " contains " << trainSet.size() << " patterns." << std::endl;
  
  // Read validation set.
  read_dataset(validFile, nInputs, nOutputs, validSet, sigma2tolerance);
  std::cerr << "Validation set: " << validFile << " contains " << validSet.size() << " patterns." << std::endl << std::endl;
  
  // Normalise training and validation sets.
  normalise();
  
  // Get initial RMS error on validation set.
  val_best = rms_error(validSet);
  
}
  

//----------------------------------------------------------------------------------
// Read in dataset from given filename.
bool Training::read_dataset(char fileName[], int nIn, int nOut, std::vector<Pattern> & data, double sigma2tolerance) {
 
  // Open data file.
  std::ifstream test_file_stream(fileName);
    
  // Temporary buffer for reading into.
  std::string buffer;
  
  // True if spectroscopic evaluation values are provided.
  bool specs;
    
  // Data.
  std::vector<double> inputs, errors, trues, sigma2s;
    
  // Reset number of records count.
  int nObjects = 0;
    
  do {
      
    // Get next non-blank line.
    do {
      getline(test_file_stream, buffer);
    } while ((buffer.find_first_not_of(' ') == std::string::npos) && (!test_file_stream.eof()));
      
    // Exit if we've reached end of file.
    if (test_file_stream.eof())
      break;
      
    // Extract data from string.
    if (!annz_util::parse_input_data(buffer, nIn, nOut, specs, inputs, errors, trues, sigma2s, sigma2tolerance))
      return false;
      
    // Check for consistency of specs.
    if (!specs) {
      std::cerr << "ERROR: Bad data file format - spectroscopic redshifts missing?" << std::endl;
      return false;
    }
      
    // Add new pattern to data vector.
    // Note errors are not used in training.
    Pattern temp;
    temp.inputs = inputs;
    temp.trues = trues;
    temp.sigma2s = sigma2s;
    data.push_back(temp);
    nObjects++;
  
  } while (true);
    
  return true;
}


//----------------------------------------------------------------
// Calculate normalisations, and apply to train and valid data sets.
void Training::normalise() {
    
  std::vector<double> sumIn(nInputs, 0.0), sumsqIn(nInputs, 0.0);
  std::vector<double> sumOut(nOutputs, 0.0), sumsqOut(nOutputs, 0.0);
  
  double tmp;
  
  //----------------------------------------------------------------------------------
  // Compute sums in preparation for calculating mean and variance of each variable.

  // Sum over training set.
  for (int pat = 0; pat < trainSet.size(); ++pat) {
    for (int inp = 0; inp < nInputs; ++inp) {
      tmp = trainSet[pat].inputs[inp];
      sumIn[inp] += tmp;
      sumsqIn[inp] += tmp * tmp;
    }
    for (int oup = 0; oup < nOutputs; ++oup) {
      tmp = trainSet[pat].trues[oup];
      sumOut[oup] += tmp;
      sumsqOut[oup] += tmp * tmp;
    }
  }

  // Sum over validation set.
  for (int pat = 0; pat < validSet.size(); ++pat) {
    for (int inp = 0; inp < nInputs; ++inp) {
      tmp = validSet[pat].inputs[inp];
      sumIn[inp] += tmp;
      sumsqIn[inp] += tmp * tmp;
    }
    for (int oup = 0; oup < nOutputs; ++oup) {
      tmp = validSet[pat].trues[oup];
      sumOut[oup] += tmp;
      sumsqOut[oup] += tmp * tmp;
    }
  }
    
  //----------------------------------------------------------------------------------------
  // Calculate means and standard deviations, and store in member variables of Train class.
  int nPatt = trainSet.size() + validSet.size();

  // Clear and assign normalisation vectors.
  meanIn.assign(nInputs, 0.0); sigmaIn.assign(nInputs, 0.0); meanOut.assign(nOutputs, 0.0); sigmaOut.assign(nOutputs, 0.0);
  
  for (int inp = 0; inp < nInputs; ++inp) {
    meanIn[inp] = (sumIn[inp]/nPatt);
    sigmaIn[inp] = sqrt(nPatt/(nPatt - 1.0) * (sumsqIn[inp]/nPatt - meanIn[inp]*meanIn[inp]));
  } 
  for (int oup = 0; oup < nOutputs; ++oup) {
    meanOut[oup] = sumOut[oup]/nPatt;
    sigmaOut[oup] = sqrt(nPatt/(nPatt - 1.0) * (sumsqOut[oup]/nPatt - meanOut[oup]*meanOut[oup]));
  }

  //for (int oup = 0; oup < nOutputs; ++oup) {
  //  std::cerr << "normalizing output #" << oup+1 << " by dividing by " << sigmaOut[oup] << std::endl;
  //}

  //----------------------------------------------------------
  // Apply normalisation to data.
  
  // Normalise training set.
  for (int pat = 0; pat < trainSet.size(); ++pat) {
    for (int inp = 0; inp < nInputs; ++inp) {
      trainSet[pat].inputs[inp] = (trainSet[pat].inputs[inp] - meanIn[inp])/sigmaIn[inp];
    }
    for (int oup = 0; oup < nOutputs; ++oup) {
      trainSet[pat].trues[oup] = (trainSet[pat].trues[oup] - meanOut[oup])/sigmaOut[oup];
      trainSet[pat].sigma2s[oup] = sigmaOut[oup]*sigmaOut[oup];
    }
  }
  
  // Normalise validation set.
  for (int pat = 0; pat < validSet.size(); ++pat) {
    for (int inp = 0; inp < nInputs; ++inp) {
      validSet[pat].inputs[inp] = (validSet[pat].inputs[inp] - meanIn[inp])/sigmaIn[inp];
    }
    for (int oup = 0; oup < nOutputs; ++oup) {
      validSet[pat].trues[oup] = (validSet[pat].trues[oup] - meanOut[oup])/sigmaOut[oup];
      validSet[pat].sigma2s[oup] /= sigmaOut[oup]*sigmaOut[oup];
    }
  }

}


//----------------------------------------------------------------------
// Write architecture, normalisation and weight values to file.
void Training::saveNetState(char file_name[]) {

  std::ofstream out_file(file_name);
  out_file << "# Generated by annz_train" << std::endl;
  out_file << "# Network architecture" << std::endl;
  
  myNetwork.dump_arch(out_file);

  out_file << std::endl;

  out_file << "# Normalisation table" << std::endl;

  for (int i = 0; i < nInputs; ++i)
    out_file << "INORM " << i << " " << meanIn[i] << " " << sigmaIn[i] << std::endl;

  for (int i = 0; i < nOutputs; ++i)
    out_file << "ONORM " << i << " " << meanOut[i] << " " << sigmaOut[i] << std::endl;
  
  out_file << std::endl;
  myNetwork.setWeights(best_wts);
  myNetwork.dump_weights(out_file); 
  out_file.close();
}
