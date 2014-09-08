#ifndef TEST_H
#define TEST_H

#include "test_network.h"

//-------------------------------------------------------------------
// Class encapsulating the committee concept - a group of networks
// with the same interface, allowing averaging of results.
//-------------------------------------------------------------------
class Committee {
  
  // Our networks.
  std::vector<TestNetwork*> nets;

public:
  
  // Results from the committee are returned in the form of one of these.
  struct Result {
    std::vector<double> outputs, errors;
  };
  
  // Committee structure.
  int n_nets;
  int nInputs, nOutputs;
  
  // Constructor.
  Committee(int n_nets, char** net_files);

  // Present given input vector to the committee - returns results in 'result'
  void pass(const std::vector<double> & inputs, 
	    const std::vector<double> & errors, 
	    Result & result);
  
};

//-----------------------------------------------------------------------
// Class controlling the testing process.
//------------------------------------------------------------------------
class Testing {

  Committee myCommittee;
  int n_nets;
  int nInputs, nOutputs;

public:

  // Constructor.
  Testing(int n_nets, char** net_files);

  // Present given testing set to committee.
  bool test(char test_file[], char out_file[]);  
  
};

#endif // TEST_H
