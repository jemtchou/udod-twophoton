#ifndef _VMDc_h
#define _VMDc_h

#include "Hadronic.hh"

namespace UDOD{
  class VMDc: public Hadronic
  {
  private:
    double mRhoSq;
    double mOmegaSq;
    double mPhiSq;
    double m0sq;
    double rRho;
    double rOmega;
    double rPhi;
    double rc;
    double ksi;

    int idchan;
    /*
      bool b_leptonic;
      bool b_hadronic;
      bool b_hadronicBW;
    */
  
  public:
    VMDc()
    {
     idchan =13;
      /*
	b_leptonic = false;
	b_hadronic = true;
	b_hadronicBW = false;
      */
    
      mRhoSq = 0.775*0.775;   //GeV^2
      mOmegaSq = 0.782*0.782; //GeV^2
      mPhiSq = 1.020*1.020;   //GeV^2
      m0sq = 1.96;            // ???, GeV^2
      rRho = 0.65;
      rOmega = 0.08;
      rPhi = 0.05;
      rc = 1-(rRho+rOmega+rPhi);
      ksi = 1./4;
    };    
    ~VMDc(){};
 
    double hT(double Q2);
    double hS(double Q2);
  
  
  };

}

#endif
