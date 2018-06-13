#ifndef _SigmaUdod_leptonic_h
#define _SigmaUdod_leptonic_h

#include <iostream>
#include "UdodChannelBase.h"

// gamma gamma -> e+e-(mu+mu-)
namespace UDOD{

  class Leptonic :  public UdodChannelBase
  {
  private:
    double m_m2;
  public:
    //m2 = me2 or mmu2;
    Leptonic(double m2) {
      idchan =1;
      m_m2 = m2;
    };    
    ~Leptonic(){};

    double x(double W2,double t1,double t2);
    double q1muq2mu(double W2,double t1,double t2);
    double Delta_t(double W2,double t1,double t2);
    double L(double t1,double t2,double W2);
    double T(double W2,double t1,double t2);
    
    //-------------sigmaAB--------------
    double sigmaTT(double W2,double t1,double t2);
    double sigmaTS(double W2,double t1,double t2);
    double sigmaST(double W2,double t1,double t2);
    double sigmaSS(double W2,double t1,double t2);
    double tauTT(double W2,double t1,double t2);
    double tauTS(double W2,double t1,double t2);
    
    double Sigma(double s,double s1,double s2,double t1,double t2,double W2,double phi,double theta);

    double SigmaTheorEE(double s,double s1,double s2,double t1,double t2,double W2,double phi,double theta);
    double SigmaTheorMuMu(double s,double s1,double s2,double t1,double t2,double W2,double phi,double theta);
    double SigmaTheorll(double s,double s1,double s2,double t1,double t2,double W2,double phi,double theta);
    double SigmaTheorEE1(double s,double s1,double s2,double t1,double t2,double W2,double phi,double theta);
  };
  
}

#endif
