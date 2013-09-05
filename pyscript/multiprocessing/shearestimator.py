import sys
import numpy as np
sys.path.append('/Users/reiko/2code/python/mylib')

class ShearEstimator:

    """ shearestimator.py
    Adopted from garyb, ShearEstimator.h.
    
    Class to accumulate individual shape measures and spit out
    the estimated shear of the population.
    Includes non-linear effects, weighting with arbitrary function,
    but not yet anything with noise involved.
    Isotropic parent population assumed.
    """

    def __init__(self):
        self.sume1 = self.sume2 = self.sumw = self.sumr1 = self.sumr3 = 0. 
        self.sumvar1 = self.sumvar2 = self.sumwesq = 0.
        self.N = 0
        self.Delta = 1.
        return
        
    def clear(self):
        self.sume1 = self.sume2 = self.sumw = self.sumr1 = self.sumr3 = 0. 
        self.sumvar1 = self.sumvar2 = self.sumwesq = 0.
        self.N = 0
        return

    def add(self, eta1, eta2, egFix=0., sig1=0., sig2=0.):
        """Add a new measurement to the ensemble.
        Note that uncertainties are assumed to be on eta1 & eta2.
        [Not yet considering correlation between them, which should avg
        to zero anyway.]
        """
        eta = np.sqrt(eta1*eta1 + eta2*eta2)
        w = 1.
        wp = wpp = wppp = 0.
        self.sume1 += w*eta1
        self.sume2 += w*eta2
        self.sumw += w
        self.sumvar1 += w*w*sig1*sig1
        self.sumvar2 += w*w*sig2*sig2
        self.sumwesq += w*w*eta*eta
        if eta == 0.:
            Y = 1.
        else:
            Y = eta / np.tanh(eta)
        # Guessing how egFix and w will interact
        # Also only keeping 2 terms of 1/egFix so noise won't explode
        # ??      sumr1 += 0.5*( (wp*eta+w)*(1+egFix+egFix*egFix) + Y*w)
        self.sumr1 += 0.5*( (wp*eta+w)/(1-egFix) + Y*w)
        if eta < 1e-3: 
            self.sumr3 += (1./16.)*(4*w/3. + 5*wpp)
        else:
            self.sumr3 += (1./16.) * ( w * (2 - Y + Y*Y*(Y-1)/(eta*eta)) + \
                                       eta*wp*(2 + (4*Y-Y*Y)/(eta*eta)) + \
                                       wpp * (2*Y+3) + \
                                       eta*wppp)
        self.N += 1
        return

    def addEnsemble(self, rhs, isLine=True):
        """Add an ensemble to an existing ensemble.
        """
        if isLine:
            rhs = rhs.split()
            if len(rhs) != 9:
                print "ShearEstimator::addEnsemble(): string of incorrect format"
                sys.exit(1)
            self.sume1 += float(rhs[0])
            self.sume2 += float(rhs[1])
            self.sumw += float(rhs[2])
            self.sumr1 += float(rhs[3])
            self.sumr3 += float(rhs[4])
            self.sumvar1 += float(rhs[5])
            self.sumvar2 += float(rhs[6])
            self.sumwesq += float(rhs[7])
            self.N += int(rhs[8])
        else:
            self.sume1 += rhs.sume1
            self.sume2 += rhs.sume2
            self.sumw += rhs.sumw
            self.sumvar1 += rhs.sumvar1
            self.sumvar2 += rhs.sumvar2
            self.sumwesq += rhs.sumwesq
            self.sumr1 += rhs.sumr1
            self.sumr3 += rhs.sumr3
            self.N += rhs.N
        return
      
    def getString(self):
        retval = '# sume1 sume2 sumw sumr1 sumr3 sumvar1 sumvar2 sumwesq n (ShearEstimator.h)\n'
        retval = '%f %f %f %f %f %f%f %f %d\n' % (self.sume1, self.sume2, self.sumw, 
                                                  self.sumr1, self.sumr3, self.sumvar1, self.sumvar2,
                                                  self.sumwesq, self.N)
        return retval

    def sigmaE(self, shapeNoise=True):
        """Return uncertainty on response-adjusted mean distortion
        Error estimate is only valid if all input shears had
        a noise-induced uncertainty specified.
        """
        if self.sumw == 0.:
            sig1 = sig2 = 0.
            return sig1, sig2

        var1 = self.sumvar1
        var2 = self.sumvar2
        if shapeNoise:
            var1 += 0.5 * self.sumwesq
            var2 += 0.5 * self.sumwesq
        sig1 = np.sqrt(var1 / (self.sumr1*self.sumr1))
        sig2 = np.sqrt(var2 / (self.sumr1*self.sumr1))
        return sig1, sig2

    def response(self):
        """Return expected linear & cubic responses of weighted mean
        """
        resp1 = self.sumr1 / self.N;
        resp3 = self.sumr3 / self.N;
        return resp1, resp3

    def shear(self):
        """Return estimated shear that produced the population
        *** note: returns distortion, not eta!
        """
        if (self.N==0 or self.sumw==0.):
            return (0., 0.)
        # Response-adjusted mean eta values:
        m1 = self.sume1 / self.sumr1
        m2 = self.sume2 / self.sumr1
        # Need to adjust the amplitude for cubic response
        M = np.hypot(m1,m2)
        if (M==0.): 
            return (0.,0.)
        G = M + pow(M,3.)*self.sumr3/self.sumr1
        G = M + pow(G,3.)*self.sumr3/self.sumr1
        G = M + pow(G,3.)*self.sumr3/self.sumr1
        G *= 0.5
        # Now G is amplitude of best-estimated reduced shear
        # Rescale into a distortion:
        m1 *= (2*G) / (1+G*G) / M;
        m2 *= (2*G) / (1+G*G) / M;
        # convert to eta
        return (m1, m2)  # This is (e1, e2) [distortion]; not (eta1, eta2)

    def meanSigEta(self):
        """Return mean and variance of weighted eta, no responsivity
        """
        if sumw == 0.:
            mean1=mean2=sig1=sig2=0.
            return mean1, mean2, sig1, sig2
        mean1 = self.sume1 / self.sumw;
        mean2 = self.sume2 / self.sumw;
        sig1 = np.sqrt((self.sumvar1/self.sumw - mean1*mean1)/self.sumw);
        sig2 = np.sqrt((self.sumvar2/self.sumw - mean2*mean2)/self.sumw);
        return mean1, mean2, sig1, sig2
        
    def setDelta(self, Delta_):
        self.Delta = Delta_
        return

    def getDelta(self):
        return self.Delta
    
    def getN(self):
        return self.N

    """
    def getW(self, eta):
        wm3 = wf(abs(eta-3*Delta/2));
        wm1 = wf(abs(eta-1*Delta/2));
        double wp1 = wf(abs(eta+1*Delta/2));
        double wp3 = wf(abs(eta+3*Delta/2));
        w = (-wp3 + 9*wp1 + 9*wm1 - wm3) / 16;
        wp = (-wp3 + 27*wp1 - 27*wm1 + wm3) / (24*Delta);
        wpp = (wp3  - wp1 - wm1 + wm3) / (2*Delta*Delta);
        wppp = (wp3 - 3*wp1 + 3*wm1 - wm3) / (Delta*Delta*Delta);
        return w, wp, wpp, wppp
    """
