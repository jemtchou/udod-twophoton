#ifndef _rho_pole0_h
#define _rho_pole0_h

#include "Hadronic.hh"

namespace UDOD{
  class rho_pole0:  public Hadronic
  {
  private:
    double mRhoSq;

    int idchan;
    /*
      bool b_leptonic;
      bool b_hadronic;
      bool b_hadronicBW;
    */
  
  public:
    rho_pole0()
    {
     idchan =12;
      /*
	b_leptonic = false;
	b_hadronic = true;
	b_hadronicBW = false;
      */ 
      mRhoSq=0.775*0.775; 
    };    
    ~rho_pole0(){};

    double hT(double Q2);

    double hS(double Q2);
  
  
  };

}

#endif
