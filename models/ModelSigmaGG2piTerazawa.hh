#ifndef _ModelSigmaGG2piTerazawa_h
#define _ModelSigmaGG2piTerazawa_h

#include "SigmaGammaGamma.hh"

namespace UDOD{
//N=5, gamma gamma -> pi+pi-, Terazawa
  class ModelSigmaGG2piTerazawa : public SigmaGammaGamma
{
private:
  double mV = 1.4; // GeV;;
  double mpi = 0.140; // GeV;
  double alpha = 1./137;
public:
  ModelSigmaGG2piTerazawa(Udod* p) : SigmaGammaGamma(p){idchan =4;
    /*
    b_leptonic = false;
    b_hadronic = true;
    b_hadronicBW = false;
    */
};    
  ~ModelSigmaGG2piTerazawa(){};

  double fpi(double x);
  double SigmaGG(double W2, double theta, double s, double s1, double s2, double t1, double t2);
};

}

#endif
