#ifndef _rho_pole_h
#define _rho_pole_h

#include "Hadronic.hh"
namespace UDOD{
  //class SigmaGammaGamma;
  
  class rho_pole: public Hadronic
  {
  private:
    double mRhoSq; 
    double ksi;
  
    int idchan;
    /*
      bool b_leptonic;
      bool b_hadronic;
      bool b_hadronicBW;
    */
  
  public:
    rho_pole()
    {
      idchan =11;
      /*
	b_leptonic = false;
	b_hadronic = true;
	b_hadronicBW = false;
      */
    
      mRhoSq=0.775*0.775; 
      ksi = 1./4;
    };    
    ~rho_pole(){};
  
    double hT(double Q2);
    double hS(double Q2);
  
  
  };

}

#endif
