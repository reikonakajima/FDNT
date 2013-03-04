// $Id: GLSimple.cpp,v 1.3 2009/11/21 04:01:52 garyb Exp $
// Split up the source code for big class

const double DEFAULT_TOLERANCE = 0.001; //  ??
const double DEFAULT_MAX_MASK_MU = 4.; // ??
const int DEFAULT_MAX_ITERATIONS = 40; // ???
const double minimumMaskSigma=3.5;  // Smallest mask size 
const double sqrtMaskFactor=1.4;    // Factor for mask size

const double MINIMUM_TRUST_RADIUS=0.03;
const double INITIAL_TRUST_RADIUS=0.4;

// Size of mask mismatch triggering flag or forcing remask
const double REMASK_THRESHOLD=0.2;	
const int MAXIMUM_REMASK=3;

// Will not search for ellipticity solution beyond this eta:
//const double MAXIMUM_ETA=2.;
const double MAXIMUM_ETA=4.;

// SV of linear solution is considered degenerate if
// an amplitude=b00 would change chisq by less than this:
const double CHISQ_DEGEN_THRESH=0.2;


#include "GLSimple1.cpp"
#include "GLSimple2.cpp"
#include "GLSimple3.cpp"

// Instantiate:

template class GLSimple<float>;
