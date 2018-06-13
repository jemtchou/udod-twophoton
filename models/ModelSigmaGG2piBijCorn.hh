#ifndef _ModelSigmaGG2piBijCorn_h
#define _ModelSigmaGG2piBijCorn_h

#include "SigmaGammaGamma.hh"
#include <complex>

namespace UDOD{
//N=2 gamma gamma -> pi+pi-, Bijnens, Cornet
class ModelSigmaGG2piBijCorn : public SigmaGammaGamma 
{
private:
  double mpi = 0.137; //GeV;
  double mK = 0.494; //GeV;
  double f = 0.134; //f=fpi, GeV;
  double L9r_p_L10r = 1.4*0.001; //(1.4 pm 0.4)* 10^(-3); 
  double alpha = 1./137;
public:
  ModelSigmaGG2piBijCorn(Udod* p) : SigmaGammaGamma(p){idchan = 5;
    /*
    b_leptonic = false;
    b_hadronic = true;
    b_hadronicBW = false;
    */
    
  };    
  ~ModelSigmaGG2piBijCorn(){};

  double Qpi(double W2);
  double QK(double W2);
  double Rea2pi(double W2);
  double Ima2pi(double W2);
  double Absa2pi(double W2);
  double beta2pi(double W2);
  double SigmaGG(double W2, double theta, double s, double s1, double s2, double t1, double t2);
  double delta(double W2);
  double SigmaGGBorn(double W2);
  
};

}

#endif
