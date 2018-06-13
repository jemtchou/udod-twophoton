#include "ModelSigmaGG2piBijCorn.hh"

#define _USE_MATH_DEFINES //for M_PI
#include <cmath>  
#include <cstdlib> // for exit(1)
#include <complex>
#include <iostream>
#include "CkTable.h"
#include "Udod.h"
#include "UdodConstants.h"

using namespace std;
using namespace UDOD;

/* 
double ModelSigmaGG2piBijCorn::Qpi(double W2)
{
  return (sqrt(W2/mpi-4)+sqrt(W2/mpi))/(sqrt(W2/mpi-4)-sqrt(W2/mpi));
}
double ModelSigmaGG2piBijCorn::QK(double W2)
{
  return (sqrt(W2/mK-4)+sqrt(W2/mK))/(sqrt(W2/mK-4)-sqrt(W2/mK));
}
*/
double ModelSigmaGG2piBijCorn::Qpi(double W2)
{
  return (sqrt(W2/(mpi*mpi)-4)+sqrt(W2/(mpi*mpi)))/(sqrt(W2/(mpi*mpi)-4)-sqrt(W2/(mpi*mpi)));
}
double ModelSigmaGG2piBijCorn::QK(double W2)
{
  return (sqrt(W2/(mK*mK)-4)+sqrt(W2/(mK*mK)))/(sqrt(W2/(mK*mK)-4)-sqrt(W2/(mK*mK)));
}

double ModelSigmaGG2piBijCorn::Rea2pi(double W2)
{
  if(Qpi(W2)>=0 && QK(W2)>=0){
    return 1+4*W2/(f*f)*(L9r_p_L10r)-
      1/(16*M_PI*M_PI*f*f)*(1.5*W2+mpi*mpi*log(Qpi(W2))*log(Qpi(W2))+0.5*mK*mK*log(QK(W2))*log(QK(W2)));
  }
  else if(Qpi(W2)<0 && QK(W2)<0){
    return 1+4*W2/(f*f)*(L9r_p_L10r)-
      1/(16*M_PI*M_PI*f*f)*(1.5*W2+mpi*mpi*(log(abs(Qpi(W2)))*log(abs(Qpi(W2)))-M_PI*M_PI)+0.5*mK*mK*(log(abs(QK(W2)))*log(abs(QK(W2)))-M_PI*M_PI));
  }
  else if(Qpi(W2)<0 && QK(W2)>=0){
    return 1+4*W2/(f*f)*(L9r_p_L10r)-
      1/(16*M_PI*M_PI*f*f)*(1.5*W2+mpi*mpi*(log(abs(Qpi(W2)))*log(abs(Qpi(W2)))-M_PI*M_PI)+0.5*mK*mK*log(abs(QK(W2)))*log(abs(QK(W2))));
  }
  else if(Qpi(W2)>=0 && QK(W2)<0){
    return 1+4*W2/(f*f)*(L9r_p_L10r)-
      1/(16*M_PI*M_PI*f*f)*(1.5*W2+mpi*mpi*log(abs(Qpi(W2)))*log(abs(Qpi(W2)))+0.5*mK*mK*(log(abs(QK(W2)))*log(abs(QK(W2)))-M_PI*M_PI));
  }
  else{
    cout<<"ModelSigmaGG2piBijCorn.cxx: Error Rea2pi(W2), W2 = "<< W2 <<endl; return 0;
  }
}

double ModelSigmaGG2piBijCorn::Ima2pi(double W2)
{
  if(Qpi(W2)>=0 && QK(W2)>=0){
    return 0;
  }
  else if(Qpi(W2)<0 && QK(W2)<0){
    return -1/(16*M_PI*M_PI*f*f)*
      2*M_PI*(mpi*mpi*log(abs(Qpi(W2)))+0.5*mK*mK*log(abs(QK(W2))));
  }
  else if(Qpi(W2)<0 && QK(W2)>=0){
    return -1/(16*M_PI*M_PI*f*f)*
      2*M_PI*mpi*mpi*log(abs(Qpi(W2)));
  }
  else if(Qpi(W2)>=0 && QK(W2)<0){
    return -1/(16*M_PI*M_PI*f*f)*
      2*M_PI*0.5*mK*mK*log(abs(QK(W2)));
  }
  else{
    cout<<"ModelSigmaGG2piBijCorn.cxx: Error Ima2pi(W2), W2 = "<< W2 <<endl; return 0;
  }
}
double ModelSigmaGG2piBijCorn::Absa2pi(double W2)
{
  return sqrt(Rea2pi(W2)*Rea2pi(W2)+Ima2pi(W2)*Ima2pi(W2));
}

double ModelSigmaGG2piBijCorn::beta2pi(double W2)
{
  return sqrt(1-4*mpi*mpi/W2);
  // beta is velocity of the pions in the c.m. frame
}

double ModelSigmaGG2piBijCorn::SigmaGG(double W2, double theta, double s, double s1, double s2, double t1, double t2)
{
  return M_PI*alpha*alpha/(2*W2)*beta2pi(W2)*
    (2*Absa2pi(W2)*Absa2pi(W2)-Rea2pi(W2)*4*beta2pi(W2)*beta2pi(W2)*sin(theta)*sin(theta)/(1-beta2pi(W2)*beta2pi(W2)*cos(theta)*cos(theta))+ 4*pow((beta2pi(W2)*sin(theta)),4)/pow(1-beta2pi(W2)*beta2pi(W2)*cos(theta)*cos(theta),2)
     );				    
}

double ModelSigmaGG2piBijCorn::delta(double W2)
{
  return (1-beta2pi(W2)*beta2pi(W2))/beta2pi(W2)*log((1+beta2pi(W2))/((1-beta2pi(W2))));				    
}

double ModelSigmaGG2piBijCorn::SigmaGGBorn(double W2)
{
  return M_PI*alpha*alpha/(2*W2)*beta2pi(W2)*
    (4*1-1*4*(2-delta(W2))+4*(2-2*delta(W2)+(1-beta2pi(W2)*beta2pi(W2))*(1+1./2*delta(W2))));
     //*15.396*pow(10,6)/(4*M_PI*M_PI);//nBn	    
}

