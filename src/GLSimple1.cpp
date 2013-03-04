// $Id: GLSimple1.cpp,v 1.2 2009/11/12 14:21:01 garyb Exp $

// Control flow and other simple routines for GLSimple

#include "GLSimple.h"

using namespace laguerre;

template <class T>
void
GLSimple<T>::invalidateData() {
  dataAreValid = false;
  fitIsValid = false;
  isSolved = false;
}

template <class T>
void 
GLSimple<T>::invalidateFit() {
  fitIsValid = false;
  isSolved = false;
}

template <class T>
void
GLSimple<T>::setBasis(Ellipse e) {
  // Does not invalidate the data unless user explicitly
  // asks to remask the data
  invalidateFit();
  basis = e;
}

template <class T>
void
GLSimple<T>::setOrder(int ord_) {
  // No automatric remasking for order change.
  invalidateFit();
  order = ord_;
  b.resize(order);
}

// Nonlinear solution options.
// Current solution is invalidated if we free
// something that used to be fixed.
template <class T>
void
GLSimple<T>::setCentering(bool f) {
  if (f && !centering) isSolved = false;
  centering = f;
}
template <class T>
void
GLSimple<T>::setDilating(bool f) {
  if (f && !dilating) isSolved = false;
  dilating = f;
}

template <class T>
void
GLSimple<T>::setShearing(bool f) {
  if (f && !shearing) isSolved = false;
  shearing = f;
}

template <class T>
bool
GLSimple<T>::reMask() {
  invalidateData();
  return acquireData();
}

template <class T>
void
GLSimple<T>::setDefaults() {
  x = y = dsig = invsig = 0;
  A = 0;
  M = 0;
  MG.clear();
  centering = shearing = dilating = true;
  flags = 0;
  invalidateData();
  tolerance = DEFAULT_TOLERANCE;
  maxMaskMu = DEFAULT_MAX_MASK_MU;
  maxIterations = DEFAULT_MAX_ITERATIONS;
  maskSigma = -1;
}

// Constructors
template <class T>
GLSimple<T>::GLSimple(list<FitExposure<T> > list_,
		   Ellipse e_, int startingOrder): 
  order(startingOrder), basis(e_), exposures(list_), b(order)
{
  setDefaults();
}

template <class T>
GLSimple<T>::GLSimple(FitExposure<T>& fe,
		   Ellipse e_, int startingOrder): 
  order(startingOrder), basis(e_),  b(order)
{
  setDefaults();
  exposures.push_back(fe);
}

template <class T>
GLSimple<T>::GLSimple(const Image<T> sci, const Image<T> wt,
		      Ellipse e_, int startingOrder, double sky_):
  order(startingOrder), basis(e_),  b(order)
{
  setDefaults();
  exposures.push_back(FitExposure<T>(sci, wt, -1, sky_, 1.));
}

template <class T>
GLSimple<T>::~GLSimple() {
  if (x) delete x;
  if (y) delete y;
  if (dsig) delete dsig;
  if (invsig) delete invsig;
  if (A) delete A;
  if (M) delete M;
  for (int i=0; i<MG.size(); i++)
    if (MG[i]) delete MG[i];
}

