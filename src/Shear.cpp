// $Id: Shear.cpp,v 1.4 2010/05/14 15:54:43 garyb Exp $
// Functions for the Shear class
#include "Shear.h"
#include <assert.h>
#include <limits>
// Convention for Shear addition is that s1 + s2 is shear by s2
// followed by shear of s1.  Note that this follows the
// notation in the methods paper (BJ02).

Shear &Shear::setE1E2(double e1in, double e2in) {
  hasMatrix = false;
  e1 = e1in;
  e2 = e2in;
  return *this;
}

Shear &Shear::setEBeta(double ein, double betain) {
  hasMatrix = false;
  e1 = ein*cos(2.*betain);
  e2 = ein*sin(2.*betain);
  return *this;
}

Shear &Shear::setEta1Eta2(double eta1in, double eta2in) {
  double scale;
  hasMatrix = false;
  // get ratio of e amplitude to eta amplitude:
  scale = sqrt(eta1in*eta1in + eta2in*eta2in);
  if (scale>0.001) scale = tanh(scale)/scale;
  else scale=1.;
  e1 = eta1in*scale;
  e2 = eta2in*scale;
  return *this;
}

Shear& Shear::setEtaBeta(double etain, double betain) {
  double e;
  hasMatrix = false;
  e = tanh(etain);
  e1 = e * cos(2.*betain);
  e2 = e * sin(2.*betain);
  return *this;
}

void Shear::getEta1Eta2(double& eta1, double& eta2) const {
  double scale;
  // get ratio of eta amplitude to e amplitude:
  scale = sqrt(e1*e1 + e2*e2);
  if (scale>0.001) scale = atanh(scale)/scale;
  else scale=1.;
  eta1 = e1*scale;
  eta2 = e2*scale;
}

void Shear::getG1G2(double& g1, double& g2) const {
  // get ratio of eta amplitude to e amplitude:
  double esq = getESq();
  double scale = (esq>1e-6) ? (1-sqrt(1-esq))/esq : 0.5;
  g1 = e1*scale;
  g2 = e2*scale;
}

void Shear::getE1E2(double& de1, double& de2) const {
  de1 = e1;
  de2 = e2;
}

Shear& Shear::setG1G2(double g1, double g2) {
  // get ratio of eta amplitude to e amplitude:
  double scale = 2./(1+g1*g1+g2*g2);
  e1 = g1*scale;
  e2 = g2*scale;
  return *this;
}

Shear& Shear::operator+=(const Shear& s2) {
  // In comparison with BJ02, the results are such that
  // s2 is applied first, then the this* (s1) shear
  double s1sq, s2sq, e1new;
  s1sq = e1*e1+e2*e2;
  if (s1sq==0.) { (*this)=s2; return *this;}

  hasMatrix = false;

  s2sq = s2.e1*s2.e1+s2.e2*s2.e2;
  assert(s1sq<=1. && s2sq<=1.);	//addition requires a realizable shear.

  double denom=1.+e1*s2.e1 + e2*s2.e2;
  if (denom==0.) {e1=e2=0.; return *this;}

  double temp = 1.-sqrt(1.-s1sq);
  e1new = e1 + s2.e1 + temp*(e1 * s2.e2 - e2 * s2.e1)*e2/s1sq;
  e2    = e2 + s2.e2 + temp*(e2 * s2.e1 - e1 * s2.e2)*e1/s1sq;
  e1 = e1new/denom;
  e2 /= denom;

  return *this;
}

Shear& Shear::operator-=(const Shear& s2) {
  // NOTE that s1 -= s2 will produce first (-s2) applied, then
  // (+s1) applied according to the local convention.
  return Shear::operator+=(-s2);
}

Shear Shear::operator+(const Shear& s2) const {
  // returns s1 + s2
  Shear out=*this;
  out += s2;   // s2 applied first, then s1
  return out;
}

Shear Shear::operator-(const Shear& s2) const {
  // returns s1 - s2
  Shear out=*this;
  out += -s2;  // -s2 applied first, then *this (s1)
  return out;
}

double
Shear::rotationWith(const Shear& rhs) const {
  double a, b, c;
  double ra, rb, rc;
  double tot11, tot21;
  getMatrix(a,b,c);
  rhs.getMatrix(ra, rb, rc);
  tot11 = a*ra + c*rc;
  tot21 = c*ra + b*rc;
  Shear sum = -(*this + rhs);
  sum.getMatrix(ra, rb, rc);
  double cc = ra * tot11 + rc * tot21;
  double ss = rc * tot11 + rb * tot21;
  return atan2(ss, cc);
}

Shear& Shear::operator*=(const double d) {
  e1 *= d; e2 *= d;
  hasMatrix = false;
  return *this;
}

Shear& Shear::operator/=(const double d) {
  hasMatrix = false;
  return Shear::operator*=(1./d);
}

Shear Shear::operator*(const double d) {
  Shear out=*this;
  out *= d;
  return out;
}

Shear Shear::operator/(const double d) {
  Shear out=*this;
  out *= 1./d;
  return out;
}

Shear operator*(const double d, const Shear& s) {
  Shear out=s;
  out *= d;
  return out;
}

void Shear::write(ostream& fout) const {
  fout << "(" << e1 << "," << e2 << ")" ;
}

ostream& 
operator<<(ostream& os, const Shear& s) {
  s.write(os);
  return os;
}

void Shear::read(istream& fin) {
  char ch;
  hasMatrix = false;
  fin >> ch >> e1 >> ch >> e2 >> ch ;
}

istream& 
operator<<(istream& is, Shear& s) {
  s.read(is);
  return is;
}

void Shear::calcMatrix() const {
  //  Matrix is defined here to be for forward point map, source plane
  // to image plane for a circular source that acquires this shape.
  // +eta in xx posn.
  if (hasMatrix) return;
  double esq= e1*e1+e2*e2;
  //Assert (esq<1.);	//Must be realizable
  // DG: sometimes floating point uncertainty causes esq==1
  // it is not clear what is the best way of dealing with this, maybe multiplying all shears by some number slighly smaller than 1
  // the Assert in any case makes the whole run abort
  // therefore I am making this very bad approximation below for now; only happens in ~1 galaxy in 1e5 or so
  if(esq>=1.) // not realizable! setting to fake and shouting out loud
  {
     cerr << "WARNING: non-realizable matrix in Shear::calcMatrix with e1=" << e1 << ", e2=" << e2 << endl;
     matrixA = 1.+0.5*e1;
     matrixB = 1.-0.5*e1;
     matrixC = +0.5*e2;
     hasMatrix = true;
     return;
  }

  
  if (esq<1e-3) {
    //Small-e approximation ok to part in 10^-6.
    matrixA = 1.+0.5*e1;
    matrixB = 1.-0.5*e1;
    matrixC = +0.5*e2;
  } else {
    double temp = sqrt(1-esq);
    double cc=sqrt(0.5*(1+1./temp));
    temp = (1-temp)/esq;
    matrixA = cc*(1+temp*e1);
    matrixB = cc*(1-temp*e1);
    matrixC = +cc*temp*e2;
  }
  hasMatrix = true;
}

DMatrix
Ellipse::getMatrix() const {
  double a, b, c;
  double scale=exp(mu);
  s.getMatrix(a,b,c);
  DMatrix m(2,2);
  m(0,0) = a*scale;
  m(1,1) = b*scale;
  m(0,1) = c*scale;
  m(1,0) = c*scale;
  return m;
}

Ellipse
Ellipse::fromMatrix(const DMatrix& m, double& rotation, bool& parityFlip) {
  Assert(m.nrows()==2 && m.ncols()==2);
  double det = m(0,0)*m(1,1) - m(0,1)*m(1,0);
  parityFlip = false;
  double scale;
  if (det < 0) {
    parityFlip = true;
    scale = -det;
  } else if (det==0.) {
    // Degenerate transformation.  Return some junk
    return Ellipse(0., 0., -std::numeric_limits<double>::max(), 0., 0.);
  } else {
    scale = det;
  }
  // Determine and remove the dilation
  double mu = 0.5*log(scale);

  // Now make m m^T matrix, which is symmetric
  // a & b are diagonal elements here
  double a = m(0,1)*m(0,1) + m(0,0)*m(0,0);
  double b = m(1,1)*m(1,1) + m(1,0)*m(1,0);
  double c = m(1,1)*m(0,1) + m(1,0)*m(0,0);
  
  double eta = acosh(MAX(1.,0.5*(a+b)/scale));
  double beta = 0.5*atan2(2.*c, a-b);
  Shear s;
  s.setEtaBeta(eta,beta);
  s.getMatrix(a,b,c);

  // Now look for the rotation
  rotation = atan2(-c*m(0,0)+a*m(1,0), b*m(0,0)-c*m(1,0));
  return Ellipse(s,mu, Position<double>(0.,0.));
}

// Ellipses share the ordering conventions:  e1 + e2 is transform
// e2 followed by transform e1.  Transform objects, not coords.
Ellipse Ellipse::operator-() const {
  Position<double> x3(-x0);
  x3 /= expmu;
  Ellipse* temp = new Ellipse(-s, -mu, s.inv(x3));
  return *temp;
}

Ellipse& Ellipse::operator+=(const Ellipse& e2) {
  Position<double> x3 = fwd(e2.getX0());
  x0 = x3;
  s += e2.getS();
  mu += e2.getMu();
  expmu = exp(mu);
  return *this;
}

Ellipse& Ellipse::operator-=(const Ellipse& e2) {
  return operator+=(-e2);
}

Ellipse Ellipse::operator+(const Ellipse& rhs) const {
  Ellipse out(*this);
  out += rhs;
  return out;
}

Ellipse Ellipse::operator-(const Ellipse& rhs) const {
  Ellipse out(*this);
  out -= rhs;
  return out;
}

void Ellipse::write(ostream& fout) const {
  s.write(fout); fout << " " << mu << " " ; x0.write(fout);
}

void Ellipse::read(istream& fin) {
  s.read(fin); fin >> mu; x0.read(fin);
}

Bounds<double>
Ellipse::range(double sig) const {
  double a,b,c;
  s.getMatrix(a,b,c);
  // ??? note that below depends on s matrix being inverse and
  // with unit determinant
  double xmax=sqrt(a*a+c*c);
  double ymax=sqrt(b*b+c*c);
  return Bounds<double>( x0.x - xmax*expmu*sig,
			 x0.x + xmax*expmu*sig,
			 x0.y - ymax*expmu*sig,
			 x0.y + ymax*expmu*sig);
}

