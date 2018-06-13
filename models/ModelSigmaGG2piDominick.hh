#ifndef _ModelSigmaGG2piDominick_h
#define _ModelSigmaGG2piDominick_h

#include "SigmaGammaGamma.hh"

namespace UDOD{
  //N=3 gamma gamma -> pi+pi-, Phys.ReV.D50,N5, Dominick et al.
  class ModelSigmaGG2piDominick : public SigmaGammaGamma 
  {
    // private:
    //double m_phi;
  public:
    // double phi;
    ModelSigmaGG2piDominick(Udod* p) : SigmaGammaGamma(p){idchan =6;
      // m_phi=phi;
    };    
    ~ModelSigmaGG2piDominick(){};
  
    //double costheta1(double theta, double s, double s1, double s2, double W2, double t1, double t2);
    double costheta1(double theta); 
    double SigmaGG(double W2, double theta, double s, double s1, double s2, double t1, double t2);
  };

}

#endif
