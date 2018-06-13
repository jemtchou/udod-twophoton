#ifndef _ModelSigmaGGpi0_h
#define _ModelSigmaGGpi0_h

#include "SigmaGammaGamma.hh"

namespace UDOD{
//N=6 gamma gamma -> pi0
class ModelSigmaGGpi0 : public SigmaGammaGamma 
{
private: 
  double alpha=1./137;
public:
  ModelSigmaGGpi0(Udod* p) : SigmaGammaGamma(p){idchan =4;
    /*
    b_leptonic = false;
    b_hadronic = true;
    b_hadronicBW = false;
    */
  };    
  ~ModelSigmaGGpi0(){};
  
  double SigmaGG(double W2, double theta, double s, double s1, double s2, double t1, double t2);
};

}

#endif
