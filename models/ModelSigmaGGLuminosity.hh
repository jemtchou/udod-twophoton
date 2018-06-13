#ifndef _ModelSigmaGGLuminosity_h
#define _ModelSigmaGGLuminosity_h

#define _USE_MATH_DEFINES //for M_PI
#include <cmath>  
#include <cstdlib> // for exit(1)
#include <complex>

#include "UdodConstants.h"

namespace UDOD{
  //N=1: Luminosity
  class ModelSigmaGGLuminosity : public SigmaGammaGamma 
  {  
  public:
    ModelSigmaGGLuminosity(Udod* p) : SigmaGammaGamma(p){idchan =9;
      /*
	b_leptonic = false;
	b_hadronic = true;
	b_hadronicBW = false;
      */
    }
    ~ModelSigmaGGLuminosity(){};

    double SigmaGG(double W2, double theta, double s, double s1, double s2, double t1, double t2)
    {
      return 1;
      //1 cm^2 = 10^24 barn;
      //1 barn = 1/0.389*10^3*1/GeV^2;
      //Budnev, p.186:
      //return alpha*alpha/(3*m_pi*m_pi);//1/GeV^2
	}
  };

}

#endif
