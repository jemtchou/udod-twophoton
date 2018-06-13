#ifndef _GVMD_h
#define _GVMD_h


#include "Hadronic.hh"

namespace UDOD{
  class GVMD: public Hadronic 
  {
  private:
    double m1sq;
    double m2sq;
    double r;
    double ksi;

    int idchan;
    /*
      bool b_leptonic;
      bool b_hadronic;
      bool b_hadronicBW;
    */
  
  public:
    GVMD() {
      idchan =55;
      /*
	b_leptonic = false;
	b_hadronic = true;
	b_hadronicBW = false;
      */
  
      m1sq = 0.54;   //GeV^2
      m2sq = 1.8; //GeV^2
      r = 3./4;
      ksi = 1./4;
    };    
    ~GVMD(){};

    double Q2;
    double P1 = 1+Q2/m1sq;
    double P2 = 1+Q2/m2sq;

    double hT(double Q2);
    double hS(double Q2);
  
  
  };

}

#endif
