// $Id: GLSimple3.cpp,v 1.6 2012/02/23 22:28:01 garyb Exp $
#include "GLSimple.h"
using namespace laguerre;


// Fitting routines

template <class T>
LVector
GLSimple<T>::linearSolution() {
  if (fitIsValid) return b;

  // Clear relevant flags first
  flags &= ~(BDegeneracy + Singularity);
  if (!acquireData()) return b;

  // Make ourselves a design matrix if we don't
  // have one of the right size:
  const int npts = dsig->size();
  const int ncoeff = PQIndex::size(order);
  DOF = npts - ncoeff;
  if (DOF<0) {
    // Not enough data:
    flags |= (BDegeneracy + Singularity);
    return b;
  }
  if (!A) {
    A = new DMatrix(npts, ncoeff);
  } else if (A->nrows()!=npts || A->ncols()!=ncoeff) {
    delete A;
    A = new DMatrix(npts, ncoeff);
  }

  // Map x and y coords into current basis:
  DVector x0 = *x; x0.addToAll(-basis.getX0().x);
  DVector y0 = *y; y0.addToAll(-basis.getX0().y);
  DMatrix mm = basis.getMatrix().inverse();
  DVector xunit = mm(0,0) * x0 + mm(0,1) * y0;
  DVector yunit = mm(1,0) * x0 + mm(1,1) * y0;
  double sigma = exp(basis.getMu());

  // Create design matrix and solve with SVD:
  LVector::design(*A, xunit, yunit, *invsig,
		  order, sigma);
  A->divideUsing(tmv::SV);
  A->saveDiv();
  try {
    A->resetDiv();
    b.rVector() = (*dsig) / (*A);
  } catch (tmv::Singular) {
    flags |= (Singularity + BDegeneracy);
    return b;
  }
  // DONE!!

  // Some heuristics to assign the degeneracy flag:
  // First, calculate the value and uncertainty that the first
  // coefficient would have (Gauss-weighted flux) if it were
  // the only thing being fit:
  var00 = 1./ ( A->col(0) * A->col(0) );
  f00 = (*dsig * A->col(0) ) * var00;

  // Then if there are b coeff's the size of f00 would not cause
  // a chisq change of CHISQ_DEGEN_THRESH, then we set the flag:
  if ( f00 * A->svd().getS()(ncoeff-1) * f00 < CHISQ_DEGEN_THRESH)
    flags |= BDegeneracy;

  // Calculate chisq of this fit:
  DVector r = *dsig * A->svd().getU();
  chisq = *dsig * *dsig - r*r;	// ??? Modify if SVD solution set some SV's to zero
  fitIsValid = true;
  return b;
}

template <class T>
void
GLSimple<T>::b00(double& value, double& var) const {
  if (fitIsValid) {
    value = f00;
    var = var00;
  } else {
    value = var = -1.;
  }
}

template <class T>
tmv::ConstVectorView<double>
GLSimple<T>::fisherSV_S() const {
  if (!A || !A->divIsSet()) 
    throw GLError("Request for Fisher matrix S that is not calculated");
  return A->svd().getS().diag();
}

template <class T>
tmv::ConstMatrixView<double>
GLSimple<T>::fisherSV_V() const {
  if (!A || !A->divIsSet()) 
    throw GLError("Request for Fisher matrix V that is not calculated");
  return A->svd().getVt();
}

template <class T>
tmv::SymMatrix<double> 
GLSimple<T>::covar() {
  if (!A || !A->divIsSet()) 
    throw GLError("Request for covar() that is not calculated");
  DMatrix tmp(A->ncols(), A->ncols());
  A->makeInverseATA(tmp);
  return tmv::SymMatrix<double>(tmp);
}

template <class T>
vector<LVector::GType>
GLSimple<T>::buildTestIndices() {
  vector<LVector::GType> v;
  if (centering) {
    v.push_back(LVector::iX);
    v.push_back(LVector::iY);
  }
  if (dilating) {
    v.push_back(LVector::iMu);
  }
  if (shearing) {
    v.push_back(LVector::iE1);
    v.push_back(LVector::iE2);
  }
  return v;
}

template <class T>
DMatrix
GLSimple<T>::dMdE(double xyscale) {
  // Need to have a successful fit to work
  if (!fitIsValid)
    throw GLError("Request for dMdE without valid current fit");

  // See if M and MG are current
  updateMG();

  // Structure the matrix elements depending on what is currently varying
  vector<LVector::GType> mIndices = buildTestIndices();
  vector<LVector::GType> eIndices = mIndices;
  if (shearing) {
    eIndices.push_back(LVector::iRot); // Need to include some rotation too
  }

  // Matrix dMde is derivs of roundness tests w.r.t. differential transform
  // from the current basis, including rotation.
  DMatrix dMde(mIndices.size(), eIndices.size());

  // Contract MG with current b vector to get dM/de
  for (int ie=0; ie<eIndices.size(); ie++) {
    int iTrans = eIndices[ie];
    DVector dM = (*MG[iTrans]) * b.rVector();
    for (int im=0; im<mIndices.size(); im++)
      dMde(im, ie) = dM[ mIndices[im]];
  }

  // Matrix dedE is derivs of differential transforms with respect to
  // the values of the basis ellipse.
  DMatrix dedE(eIndices.size(), mIndices.size(), 0.);

  // Prepare some things:
  double eta1, eta2;
  basis.getS().getEta1Eta2(eta1, eta2);
  double eta = basis.getS().getEta();
  double Y1 = eta>0. ? sinh(eta) / eta : 1.;
  double Y2 = eta>0. ? pow( sinh(eta/2.) / eta, 2.) : 0.25;

  double sa, sb, sc;
  basis.getS().getMatrix(sb,sa,sc);
  sa *= exp(-basis.getMu()) * xyscale;
  sb *= exp(-basis.getMu()) * xyscale;
  sc *= -exp(-basis.getMu()) * xyscale;

  for (int ie=0; ie<eIndices.size(); ie++) {
    LVector::GType iDiff = eIndices[ie];
    for (int im=0; im<mIndices.size(); im++) {
      LVector::GType iGlob = mIndices[im];

      // Double-dispatch
      if (iGlob==LVector::iMu) {
	// Dilation has one component on diagonal
	if (iDiff==LVector::iMu)
	  dedE(ie, im) = 1.;
      } else if (iGlob==LVector::iE1) {
	// Each shear component influence 2 shears + rotation:
	if (iDiff==LVector::iE1)
	  dedE(ie, im) = eta>0. ? Y1 + (1-Y1)*eta1*eta1/(eta*eta) : 1.;
	if (iDiff==LVector::iE2)
	  dedE(ie, im) = eta>0. ? (1-Y1)*eta1*eta2/(eta*eta) : 0.;
	if (iDiff==LVector::iRot)
	  dedE(ie, im) = -Y2*eta2;
      } else if (iGlob==LVector::iE2) {
	if (iDiff==LVector::iE1)
	  dedE(ie, im) = eta>0. ? (1-Y1)*eta1*eta2/(eta*eta) : 0.;
	if (iDiff==LVector::iE2)
	  dedE(ie, im) = eta>0. ? Y1 + (1-Y1)*eta2*eta2/(eta*eta) : 1.;
	if (iDiff==LVector::iRot)
	  dedE(ie, im) = Y2*eta1;
      } else if (iGlob==LVector::iX) {
	if (iDiff==LVector::iX) 
	  dedE(ie, im) = sa;
	if (iDiff==LVector::iY) 
	  dedE(ie, im) = sc;
      } else if (iGlob==LVector::iY) {
	if (iDiff==LVector::iX) 
	  dedE(ie, im) = sc;
	if (iDiff==LVector::iY) 
	  dedE(ie, im) = sb;
      } else
	throw GLError("dMdE() encountered unknown transformation code");
    }
  }

  return -(dMde * dedE);
}

template <class T>
void
GLSimple<T>::updateM() {
  const int nTests=5;
  const int ncoeff = PQIndex::size(order);
  // We're ok if we already have M of correct order:
  if (M && M->ncols()==ncoeff) return;

  if (M) delete M;
  M = new DMatrix(nTests, ncoeff, 0.);
  // Simple tests:
  PQIndex pq(1,0);
  (*M)(LVector::iX, pq.rIndex()) = 1.;
  (*M)(LVector::iY, pq.rIndex()+1) = 1.;
  pq.setPQ(1,1);
  (*M)(LVector::iMu, pq.rIndex()) = 1.;
  pq.setPQ(2,0);
  (*M)(LVector::iE1, pq.rIndex()) = 1.;
  (*M)(LVector::iE2, pq.rIndex()+1) = 1.;
}

template <class T>
void
GLSimple<T>::updateMG() {
  updateM();
  const int ncoeff = PQIndex::size(order);
  // We're ok if we already have correct order:
  if (MG.size()==LVector::nGen && MG[0]->ncols()==ncoeff) return;

  // Clear out old ones:
  for (int i=0; i<MG.size(); i++)
    delete MG[i];
  MG.clear();

  // Start making new guys:
  for (LVector::GType i = LVector::GType(0); 
       i<LVector::nGen;
       ++i) 
    MG.push_back(new DMatrix( *M * LVector::Generator(i, order, order)));

  return;
}

template <class T>
DVector
GLSimple<T>::Mtests() {
  if (!fitIsValid)
    throw GLError("Request for Mtests() without valid current fit");

  updateM();

  DVector allM = *M * b.rVector();
  // output elements that are varying:
  vector<LVector::GType> v = buildTestIndices();
  DVector outM(v.size());
  for (int i=0; i<v.size(); i++)
    outM[i] = allM[v[i]];
  return DVector(outM);
}

template <class T>
bool
GLSimple<T>::solve() {

  // Solve with Dogleg trust region method, using current order
  // and current basis as starting point.  Use analytic
  // derivatives and see what happens.

  // Re-mask at start of the fit and if the basis changes by more
  // than REMASK_THRESHOLD in any direction, up to MAXIMUM_REMASK
  // times.  

  // Exceed MAXIMUM_ETA halts iteration and sets flag for this.
  // Exceed maxIterations also halts and sets convergence flag, as
  // does a matrix singularity or a trust region that gets too small.

  // Final solution with mismatched mask sets flag too.

  // A successful step that is below threshold size signals success.
  // ?? Include gainFactor requirement??

  if (isSolved) return true;

  // Clear relevant flags first
  flags &= ~(BDegeneracy + Singularity); // ??? more...

  // x & y shifts will be in units of xyscale so that all 5
  // ellipse components change by order unity.
  double xyscale = exp(basis.getMu());

  vector<LVector::GType> whatFit=buildTestIndices();
  const int nFit = whatFit.size();

  bool restart = true;
  int nRemasks=0;

  double eta1,eta2, x0, y0, mu;
  double trustRadius = INITIAL_TRUST_RADIUS;

  DVector bestE(nFit);
  double bestMsq;
  DVector f(nFit);
  DMatrix J(nFit, nFit);

  for (int istep=0; istep<maxIterations; istep++) {
    if (restart) {
      trustRadius = INITIAL_TRUST_RADIUS;
      if (!reMask()) return false;
      nRemasks++;
      dbg << "Solving at basis (before initial linearSolution) " << basis << endl;
      linearSolution();
      if (!fitIsValid) return false;
      // Set up the data vector from current basis
      basis.getS().getEta1Eta2(eta1,eta2);
      for (int i=0; i<bestE.size(); i++) {
	if (whatFit[i]==LVector::iX) bestE[i] = basis.getX0().x / xyscale;
	else if (whatFit[i]==LVector::iY) bestE[i] = basis.getX0().y / xyscale;
	else if (whatFit[i]==LVector::iMu) bestE[i] = basis.getMu();
	else if (whatFit[i]==LVector::iE1) bestE[i]=eta1;
	else if (whatFit[i]==LVector::iE2) bestE[i]=eta2;
      }
      // ??? set bestE and bestMsq
      f=Mtests();
      J=dMdE(xyscale);
      bestMsq = f*f;
      restart = false;
      dbg << "Start Msq " << bestMsq << " at " << bestE << endl;
    }

    // Get the value and derivatives of tests at this point:
    dbg << "*******" << endl;

    // Choose step: first is Newton-Raphson
    DVector dE(nFit);
    DVector nr(nFit);
    try {
      nr = -f/J;
    } catch (tmv::Singular) {
      flags |= (Singularity + DidNotConverge);
      return false;
    }
    dbg << "nr: " << nr << endl;
    if (nr*nr < trustRadius * trustRadius) {
      dE = nr;
    } else {
      // Calculate steepest descent vector:
      DVector grad = -J.transpose() * f;
      DVector tmp = J*grad;
      grad *= (grad*grad) / (tmp*tmp);
      dbg << "grad: " << grad << endl;
      if (grad*grad > trustRadius*trustRadius)
	// use shortened steepest descent
	dE = grad*trustRadius / sqrt(grad*grad);
      else {
	// Find where dogleg crosses trustRadius
	DVector vb = nr - grad;
	double aa = grad*grad;
	double ab = grad*vb;
	double bb = vb*vb;
	double alpha = (-ab+sqrt(ab*ab+bb*(trustRadius*trustRadius-aa)))/bb;
	Assert( alpha>=0. && alpha<=1.);
	dE = grad + alpha*vb;
      }
    }
    dbg << "dE: " << dE <<endl;
    // Set basis to the suggested trial:
    DVector tryE = bestE + dE;
    for (int i=0; i<tryE.size(); i++) {
      if (whatFit[i]==LVector::iX) x0=xyscale*tryE[i];
      else if (whatFit[i]==LVector::iY) y0=xyscale*tryE[i];
      else if (whatFit[i]==LVector::iMu) mu=tryE[i];
      else if (whatFit[i]==LVector::iE1) eta1=tryE[i];
      else if (whatFit[i]==LVector::iE2) eta2=tryE[i];
    }
    if (centering) basis.setX0(Position<double>(x0,y0));
    if (dilating) basis.setMu(mu);
    if (shearing) {
      Shear S; S.setEta1Eta2(eta1,eta2);
      basis.setS(S);
    }
    invalidateFit();

    // Quit if trial is out of bounds
    if (basis.getS().getEta() > MAXIMUM_ETA) {
      flags |= TooElliptical;
      return false;
    }

    // Restart for mask drift
    if (nRemasks < MAXIMUM_REMASK) {
      double smallest, largest;
      maskRatios(smallest, largest);
      if (largest > 1 + REMASK_THRESHOLD
	  || smallest < 1 - REMASK_THRESHOLD) {
	nRemasks++;
	// Start over with new mask at current basis
	restart = true;
	continue;
      }
    }

    // Get expected reduction
    DVector expectedDM = J*dE;
    double expectedDMsq = -2*(f*expectedDM) - expectedDM*expectedDM;

    // Make new measurement:
    dbg << "Solving at basis " << basis << endl;
    linearSolution();
    if (!fitIsValid) return false;

    // Get gain factor
    double Msq;
    DVector Mnow=Mtests();
    Msq = Mnow*Mnow;
    double gainFactor = (bestMsq - Msq) / expectedDMsq;

    dbg << " iter " << istep
	<< " at " << tryE << endl;
    dbg << " bestMsq " << bestMsq << " expected "
	<< expectedDMsq << endl;
    dbg << "   Msq " << Msq << " gainFactor " << gainFactor 
	<< "  trust " << trustRadius << endl;
    // Was there improvement??
    if (gainFactor > 0.) {
      // Mark new best point.
      bestE = tryE;
      bestMsq = Msq;
      // And save function and Hessian here
      f = Mnow;
      J = dMdE(xyscale);
      // If the movement was smaller than tolerance, we have succeeded!
      if (dE*dE < tolerance*tolerance) {
	// SUCCESS: FINISH AND EXIT
	isSolved = true;
	// Check whether we finished with a good mask fit
	double smallest, largest;
	maskRatios(smallest, largest);
	if (largest > 1 + REMASK_THRESHOLD
	    || smallest < 1 - REMASK_THRESHOLD) 
	  flags |= MaskMismatch;
	return true;
      }
    } 

    // Adjust the trustRadius
    dbg << "gainFactor: " << gainFactor << endl;
    if (gainFactor < 0.25) {
      trustRadius *= 0.5;
    } else if (gainFactor > 0.75) {
      trustRadius = MAX(trustRadius, 3.*sqrt(dE*dE));
      trustRadius = MIN(trustRadius, INITIAL_TRUST_RADIUS);
    }

    // Abort for trust radius too small
    if (trustRadius < MINIMUM_TRUST_RADIUS) {
      dbg << "trust radius too small, " << trustRadius << endl;
      flags |= DidNotConverge;
      return false;
    }
  }
  // Here for too many steps?
  flags |= DidNotConverge;
  dbg << "too many steps" << endl;
  return false;
}

