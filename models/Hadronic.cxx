
#define _USE_MATH_DEFINES //for M_PI
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <numeric>
#include <cstdlib> // for exit(1)

#include "Hadronic.hh"
#include "CkTable.h"
#include "Udod.h"
#include "UdodConstants.h"

using namespace std;
using namespace UDOD;

double Hadronic::Sigma(double s,double s1,double s2,double t1,double t2,double W2,double phi,double theta)
{
   double u1 = s2-me2-t1;
   double u2 = s1-me2-t2;
   double nu = 1./2*(W2-t1-t2); 
   double K = 1/sqrt(W2)*sqrt(nu*nu-t1*t2);
   //1: 
   double rho1PP = 1./2*((u2-nu)*(u2-nu)/(K*K*W2)+1+4*me2/t1);
   double rho100 = (u2-nu)*(u2-nu)/(K*K*W2)-1;
   double absrho1PM = rho1PP-1;
   double absrho1P0 = (u2-nu)/(K*sqrt(W2))*sqrt(rho1PP-1);
   //double absrho1PM =  abs(rho1PP-1);
   //double absrho1P0 =  abs((u2-nu)/(K*sqrt(W2))*sqrt(rho1PP-1));
   //2:
   double rho2PP = 1./2*((u1-nu)*(u1-nu)/(K*K*W2)+1+4*me2/t2);
   double rho200 = (u1-nu)*(u1-nu)/(K*K*W2)-1;
   double absrho2PM = rho2PP-1;
   double absrho2P0 =(u1-nu)/(K*sqrt(W2))*sqrt(rho2PP-1);
   //double absrho2PM =  abs(rho2PP-1);
   //double absrho2P0 =  abs((u1-nu)/(K*sqrt(W2))*sqrt(rho2PP-1));

   // phi
   double cosphi1 = cos(phi)+sin(phi)*sin(phi)*(sqrt((-t1)*(-t2))*(2*s-s1-s2))/((W2-t1-t2)*sqrt((s-s1)*(s-s2)));	
   double cos2phi1 = 2*cosphi1*cosphi1-1;
  
    return  2*rho1PP*2*rho2PP*sigmaTT(-t1,-t2,W2,theta,s,s1,s2)+
      2*rho1PP*rho200*sigmaTS(-t1,-t2,W2,theta,s,s1,s2)+
      rho100*2*rho2PP*sigmaST(-t1,-t2,W2,theta,s,s1,s2)+
      rho100*rho200*sigmaSS(-t1,-t2,W2,theta,s,s1,s2)+
      2*absrho1PM*absrho2PM*tauTT(-t1,-t2,W2,theta)*cos2phi1-
      8*absrho1P0*absrho2P0*tauTS(-t1,-t2,W2,theta)*cosphi1;
}
/*
double Hadronic::sigmaTT(double Q12, double Q22,double W2,double theta,double s, double s1, double s2){
  return hT(Q12)*hT(Q22)*psgg->SigmaGG(W2,theta,s,s1,s2,-Q12,-Q22);
}
*/

double Hadronic::sigmaTT(double Q12, double Q22,double W2,double theta,double s, double s1, double s2){
  if(pSGG==0){cout<< "No model for Gamma Gamma" << endl;}
  return hT(Q12)*hT(Q22)*pSGG->SigmaGG(W2,theta,s,s1,s2,-Q12,-Q22);
}
double Hadronic::sigmaTS(double Q12, double Q22,double W2,double theta,double s, double s1, double s2){
  if(pSGG==0){cout<< "No model for Gamma Gamma" << endl;}
  return hT(Q12)*hS(Q22)*pSGG->SigmaGG(W2,theta,s,s1,s2,-Q12,-Q22);
}
double Hadronic::sigmaST(double Q12, double Q22,double W2,double theta,double s, double s1, double s2){
  if(pSGG==0){cout<< "No model for Gamma Gamma" << endl;}
  return hS(Q12)*hT(Q22)*pSGG->SigmaGG(W2,theta,s,s1,s2,-Q12,-Q22);
}
double Hadronic::sigmaSS(double Q12, double Q22,double W2,double theta,double s, double s1, double s2){
  if(pSGG==0){cout<< "No model for Gamma Gamma" << endl;}
  return hS(Q12)*hS(Q22)*pSGG->SigmaGG(W2,theta,s,s1,s2,-Q12,-Q22);
}


double Hadronic::tauTT(double Q12, double Q22,double W2,double theta){return 0;}
double Hadronic::tauTS(double Q12, double Q22,double W2,double theta){return 0;}

// Budnev, p.186;
double Hadronic::SigmaLuminosityTheor(double s)
{

  return pow(alpha,4)/(18*M_PI*M_PI*m_pi*m_pi)*
    log(s*m_rho*m_rho/(me2*m_pi*m_pi))*
    log(s*pow(m_rho,6)/(me2*me2*me2*m_pi*m_pi))*
    log(s/(m_pi*m_pi))*log(s/(m_pi*m_pi))*
    nbarn;
}

