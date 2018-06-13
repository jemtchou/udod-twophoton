#include "ModelSigmaGG2KDominick.hh"
#include "CkTable.h"
#include "Udod.h"
#include "UdodConstants.h"

using namespace std;
using namespace UDOD;

#define _USE_MATH_DEFINES //for M_PI
#include <cmath>  

/*
double ModelSigmaGG2KDominick::costheta1(double theta, double s, double s1, double s2, double W2, double t1, double t2)
{
  return cos(theta)+sin(theta)*sin(theta)*(sqrt((-t1)*(-t2))*(2*s-s1-s2))/((W2-t1-t2)*sqrt((s-s1)*(s-s2)));			    
}
*/
double ModelSigmaGG2KDominick::costheta1(double theta)
{
  return theta;
}

// theta1 is in gamma gamma c.m. system !
double ModelSigmaGG2KDominick::SigmaGG(double W2, double theta, double s, double s1, double s2, double t1, double t2)
{
  return 4*2.2*(0.4/W2)*(0.4/W2)/(1+pow(costheta1(theta),4))*sigma2mu;
}

