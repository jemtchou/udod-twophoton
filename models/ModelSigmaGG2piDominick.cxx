#include "ModelSigmaGG2piDominick.hh"

#define _USE_MATH_DEFINES //for M_PI
#include <cmath>  
#include <cstdlib> // for exit(1)
#include <complex>

#include "CkTable.h"
#include "Udod.h"
#include "UdodConstants.h"

#include "Leptonic.hh"
//#include "UdodChannelBase.h"

using namespace std;
using namespace UDOD;

/*
double ModelSigmaGG2piDominick::costheta1(double theta, double s, double s1, double s2, double W2, double t1, double t2)
{
   return cos(theta)+sin(theta)*sin(theta)*(sqrt((-t1)*(-t2))*(2*s-s1-s2))/((W2-t1-t2)*sqrt((s-s1)*(s-s2)));			    
}
*/

double ModelSigmaGG2piDominick::costheta1(double theta)
{
  return theta;
}

// theta1 is in gamma gamma c.m. system !

double ModelSigmaGG2piDominick::SigmaGG(double W2, double theta, double s, double s1, double s2, double t1, double t2)
{
  /*
  Leptonic LeptonicMuMu(Udod* p,long double mmu2);
  double sigma2mu =1.;
  cout << "LeptonicMuMu.Sigma: " << LeptonicMuMu.Sigma(s,s1,s2,t1,t2,W2,phi,theta) << endl;
  */
  double sigma2mu = 1;
  return 4*(0.4/W2)*(0.4/W2)/(1+pow(costheta1(theta),4))*sigma2mu;
}

