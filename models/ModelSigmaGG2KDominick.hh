#ifndef _ModelSigmaGG2KDominick_h
#define _ModelSigmaGG2KDominick_h

#include "SigmaGammaGamma.hh"

namespace UDOD{
//N=4 gamma gamma -> K+K-, Phys.ReV.D50,N5, Dominick et al.
class ModelSigmaGG2KDominick : public SigmaGammaGamma 
{
private:
  double sigma2mu = 1;
public:
  ModelSigmaGG2KDominick(Udod* p) : SigmaGammaGamma(p){idchan =7;
    /*
    b_leptonic = false;
    b_hadronic = true;
    b_hadronicBW = false;
    */
  };    
  ~ModelSigmaGG2KDominick(){};
  
  //double costheta1(double theta, double s, double s1, double s2, double W2, double t1, double t2);
  double costheta1(double theta); 
  double SigmaGG(double W2, double theta, double s, double s1, double s2, double t1, double t2);
};

}

#endif
