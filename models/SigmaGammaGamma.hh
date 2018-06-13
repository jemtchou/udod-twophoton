#ifndef _SigmaGammaGamma_h
#define _SigmaGammaGamma_h

#define _USE_MATH_DEFINES //for M_PI
#include <cmath>  
#include <cstdlib> // for exit(1)
#include <complex>
#include "Udod.h"

namespace UDOD{
// SigmaGammaGamma------------------------
class SigmaGammaGamma
{
private:
  //  int NModel;
public:
  SigmaGammaGamma(Udod* p){/*NModel=0;*/};//(int NModel_){NModel=NModel_;};    
  ~SigmaGammaGamma(){};

  //virtual double SigmaGG(double W2, double theta)=0;
  virtual double SigmaGG(double W2, double theta, double s, double s1, double s2, double t1, double t2)=0;

protected:
  int idchan;
  
  //bool b_leptonic,b_hadronic,b_hadronicBW; 
  
  // virtual double SigmaGG2pi(double W2, double theta)=0;
  // virtual double SigmaGG2K(double W2, double theta)=0;
  // virtual double SigmaGG2piTerazawa(double W2, double theta)=0;

  //void SetNModel(int NModel_){NModel=NModel_;}; 
  
  //double GetCrossSection(){return 0;}
  
  /*double sigmaGG(double W2, double theta)
  {
    if(NModel==1)
      {
	return 1;
	// Luminosity
      }
    else if(NModel==2)
      {
	//gamma gamma -> pi+ pi-, Bijnens, Cornet
	return SigmaGG(W2, theta);
      }
  else if(NModel==3)
      {
	//2pi and 2K
	return SigmaGG(W2, theta);
	//return  SigmaGG2K(W2, theta);
      }
  else if(NModel==4)
      {
	//2pi Terazawa
	return SigmaGG(W2, theta);
      }
    else
      {
	// Mat
	std::cout<< "Model not defined!"<<std::endl; return 0;
      }
      }*/

};

}

#endif
