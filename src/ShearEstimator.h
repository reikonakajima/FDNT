// $Id: ShearEstimator.h,v 1.4 2010/04/26 17:09:51 garyb Exp $
// Class to accumulate individual shape measures and spit out
// the estimated shear of the population.
// Includes non-linear effects, weighting with arbitrary function,
// but not yet anything with noise involved.
// Isotropic parent population assumed.

#ifndef SHEARESTIMATOR_H
#define SHEARESTIMATOR_H

#include "Shear.h"

namespace laguerre {

  template <class W>
  class ShearEstimator {
  public:
    ShearEstimator(const W& wf_, double Delta_=0.02): wf(wf_),
						     Delta(Delta_) {clear();}
    void clear() {
      sume1=sume2=sumw=sumr1=sumr3=0.; 
      sumvar1 = sumvar2 = sumwesq = 0.;
      n=0;
    }
    ~ShearEstimator() {}
    // Add a new measurement to the ensemble.
    // Note that uncertainties are assumed to be on eta1 & eta2.
    // [Not yet considering correlation between them, which should avg
    //  to zero anyway.]
    void add(Shear S, double egFix=0., double sig1=0., double sig2=0.) {
      double eta1, eta2, eta;
      S.getEta1Eta2(eta1, eta2);
      eta = S.getEta();
      double w, wp, wpp, wppp;
      getW(eta, w, wp, wpp, wppp);
      sume1 += w*eta1;
      sume2 += w*eta2;
      sumw += w;
      sumvar1 += w*w*sig1*sig1;
      sumvar2 += w*w*sig2*sig2;
      sumwesq += w*w*eta*eta;
      double Y = (eta == 0.) ? 1. : eta / tanh(eta);
      // Guessing how egFix and w will interact
      // Also only keeping 2 terms of 1/egFix so noise won't explode
      //??      sumr1 += 0.5*( (wp*eta+w)*(1+egFix+egFix*egFix) + Y*w);
      sumr1 += 0.5*( (wp*eta+w)/(1-egFix) + Y*w);
      if (eta < 1e-3) 
	sumr3 += (1./16.)*(4*w/3. + 5*wpp);
      else
	sumr3 += (1./16.) * ( w * (2 - Y + Y*Y*(Y-1)/(eta*eta)) +
			      eta*wp*(2 + (4*Y-Y*Y)/(eta*eta)) +
			      wpp * (2*Y+3) +
			      eta*wppp);
      n++;
    }
    // Add an ensemble to an existing ensemble.
    void addEnsemble(ShearEstimator rhs){
      sume1 += rhs.sume1;
      sume2 += rhs.sume2;
      sumw += rhs.sumw;
      sumvar1 += rhs.sumvar1;
      sumvar2 += rhs.sumvar2;
      sumwesq += rhs.sumwesq;
      sumr1 += rhs.sumr1;
      sumr3 += rhs.sumr3;
      n += rhs.n;
    }
    // Return uncertainty on response-adjusted mean distortion
    // Error estimate is only valid if all input shears had
    // a noise-induced uncertainty specified.
    void sigmaE(double& sig1, double& sig2, bool shapeNoise=true) const {
      if (sumw==0.) {
	sig1=sig2=0.;
	return;
      }
      double var1 = sumvar1;
      double var2 = sumvar2;
      if (shapeNoise) {
	var1 += 0.5*sumwesq;
	var2 += 0.5*sumwesq;
      }
      sig1 = sqrt(var1 / (sumr1*sumr1));
      sig2 = sqrt(var2 / (sumr1*sumr1));
    }
    // Return expected linear & cubic responses of weighted mean
    void response(double & resp1, double& resp3) {
      resp1 = sumr1 / n;
      resp3 = sumr3 / n;
    }
    // Return estimated shear that produced the population
    operator const Shear() const {
      if (n==0 || sumw==0.) {
	return Shear(0.,0.);
      }
      // Response-adjusted mean eta values:
      double m1 = sume1 / sumr1;
      double m2 = sume2 / sumr1;
      // Need to adjust the amplitude for cubic response
      double M = hypot(m1,m2);
      if (M==0.) return Shear(0.,0.);
      double G = M + pow(M,3.)*sumr3/sumr1;
      G = M + pow(G,3.)*sumr3/sumr1;
      G = M + pow(G,3.)*sumr3/sumr1;
      G *= 0.5;
      // Now G is amplitude of best-estimated reduced shear
      // Rescale into a distortion:
      m1 *= (2*G) / (1+G*G) / M;
      m2 *= (2*G) / (1+G*G) / M;
      return Shear(m1,m2);
    }
    // Return mean and variance of weighted eta, no responsivity
    void meanSigEta(double& mean1, double& mean2, double& sig1, double& sig2) const {
      if (sumw==0.) {
	mean1=mean2=sig1=sig2=0.;
	return;
      }
      mean1 = sume1 / sumw;
      mean2 = sume2 / sumw;
      sig1 = sqrt((sumvar1/sumw - mean1*mean1)/sumw);
      sig2 = sqrt((sumvar2/sumw - mean2*mean2)/sumw);
    }
    void setDelta(double Delta_) {Delta=Delta_;}
    double getDelta() const {return Delta;}
    int getN() const {return n;}
    string getLine(bool allComment=false) const {
      ostringstream oss;
      oss << "# sume1 sume2 sumw sumr1 sumr3 sumvar1 sumvar2 sumwesq n (ShearEstimator.h)"
	  << endl;
      if (allComment) {
	oss << "# ";
      }
      oss << sume1 << " "
	  << sume2 << " "
	  << sumw << " "
	  << sumr1 << " "
	  << sumr3 << " "
	  << sumvar1 << " "
	  << sumvar2 << " "
	  << sumwesq << " "
	  << n << " "
	  << endl;
      return oss.str();
    }
  private:
    const W& wf;	// weight function
    double Delta; // Interval for getting weight derivatives
    int n;
    double sume1;
    double sume2;
    double sumw;
    double sumr1;
    double sumr3;
    double sumvar1;
    double sumvar2;
    double sumwesq;
    void getW(double eta, double& w, double& wp, double& wpp, double& wppp) const {
      double wm3 = wf(abs(eta-3*Delta/2));
      double wm1 = wf(abs(eta-1*Delta/2));
      double wp1 = wf(abs(eta+1*Delta/2));
      double wp3 = wf(abs(eta+3*Delta/2));
      w = (-wp3 + 9*wp1 + 9*wm1 - wm3) / 16;
      wp = (-wp3 + 27*wp1 - 27*wm1 + wm3) / (24*Delta);
      wpp = (wp3  - wp1 - wm1 + wm3) / (2*Delta*Delta);
      wppp = (wp3 - 3*wp1 + 3*wm1 - wm3) / (Delta*Delta*Delta);
    }
  };

  class One {
  public:
    double operator()(double eta) const {return 1.;}
  };

  class UnweightedShearEstimator: public ShearEstimator<One> {
  public:
    UnweightedShearEstimator(double Delta=0.02): localwf(new One),
						 ShearEstimator<One>(*localwf, Delta) {};
    ~UnweightedShearEstimator() {delete localwf;}
  private:
    One* localwf;
  };

} // namespace laguerre
#endif
