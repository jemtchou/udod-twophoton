#ifndef _Hadronic_h
#define _Hadronic_h

#define _USE_MATH_DEFINES //for M_PI
#include <iostream>
#include <cmath>  
#include <cstdlib> // for exit(1)

#include "SigmaGammaGamma.hh"

namespace UDOD{ 
  class Hadronic : public UdodChannelBase
  {
  private:
    SigmaGammaGamma* pSGG;

  
  public:
    Hadronic() { pSGG=0; };
      /*
      //pSGG = new ModelSigmaGGLuminosity();
      //pSGG = new ModelSigmaGG2piBijCorn();
      //pSGG = new ModelSigmaGG2piDominick();
      //pSGG = new ModelSigmaGG2KDominick();
      pSGG = new ModelSigmaGG2piTerazawa(p);
      //pSGG = new ModelSigmaGGpi0();
      */
    /*
      if(Nh==1){
	pSig = new GVMD(p,pSGG);
      }
      else if(Nh==2){
	pSig = new VMDc(p,pSGG);
      }
      else if(Nh==3){
	pSig = new rho_pole(p,pSGG);
      }
      else if(Nh==4){
	pSig = new rho_pole0(p,pSGG);
      }
      else{
	std::cout << "Error: No model " << std::endl;
	std::exit(1);
      }
      */
   
    ~Hadronic(){};
    //double SigmaHadronic(double s,double s1,double s2,double t1,double t2,double W2,double phi,double theta,double m2);

    void SetModel(SigmaGammaGamma* sgg){pSGG = sgg;}
   
    
    double Sigma(double s,double s1,double s2,double t1,double t2,double W2,double phi,double theta);
 
    double sigmaTT(double Q12, double Q22,double W2,double theta,double s, double s1, double s2);
    double sigmaTS(double Q12, double Q22,double W2,double theta,double s, double s1, double s2);
    double sigmaST(double Q12, double Q22,double W2,double theta,double s, double s1, double s2);
    double sigmaSS(double Q12, double Q22,double W2,double theta,double s, double s1, double s2);
    double tauTT(double Q12, double Q22,double W2,double theta);
    double tauTS(double Q12, double Q22,double W2,double theta);
  
    virtual double hT(double Q2)=0;
    virtual double hS(double Q2)=0;

    double SigmaLuminosityTheor(double s);

  };

}

#endif
