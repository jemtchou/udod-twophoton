#include <iostream>
#include <iomanip>
#include <vector>
#include <numeric>

#include "SigmaBW.hh"
#include "Udod.h"
#include "UdodConstants.h"

using namespace std;
using namespace UDOD;

double SigmaBW::GammaGG()
{
  //1:jp0P, 2:jp0M, 3: jp1P, 4: jp2P, 5: jp2M 
  if(m_JP==1){return 4;}
  else if(m_JP==2){return 144;}
  else if(m_JP==3){return 32;}
  else if(m_JP==4){return 192./5;}
  else if(m_JP==5){return 64;}
  else{ std::cout << "SigmaBW.cxx: Error JP " << m_JP << std::endl; return 0;}
}

double SigmaBW::nu(double W2,double t1,double t2)
{
  return 0.5*(W2-t1-t2);
}

// fAB-----------------------------
double SigmaBW::X(double W2,double Q12,double Q22)
{
  return nu(W2,-Q12,-Q22)*nu(W2,-Q12,-Q22)-Q12*Q22;//Q12=-t1, Q22=-t2;
}

double SigmaBW::kappa(double W2,double Q12,double Q22)
{
  return m_Mp*m_Mp/(2*sqrt(X(W2,Q12,Q22)));
 
}
  
// There was M2, but really it is Mp*Mp, in kappa also.
// f_TT, J=0-, 0+, 1+, 2+, 2-
double SigmaBW::fTT_0M(double W2,double Q12,double Q22)
{
  return kappa(W2,Q12,Q22)*X(W2,Q12,Q22)/(nu(W2,-Q12,-Q22)*nu(W2,-Q12,-Q22));
}

double SigmaBW::fTT_0P(double W2,double Q12,double Q22)
{
  return kappa(W2,Q12,Q22)*pow((X(W2,Q12,Q22)+nu(W2,-Q12,-Q22)*m_Mp*m_Mp)/(3*nu(W2,-Q12,-Q22)*nu(W2,-Q12,-Q22)),2);
}

double SigmaBW::fTT_1P(double W2,double Q12,double Q22)
{
  return  kappa(W2,Q12,Q22)*pow((Q22-Q12)/(2*nu(W2,-Q12,-Q22)),2);
}
  
double SigmaBW::fTT_2P(double W2,double Q12,double Q22)
{
  return  kappa(W2,Q12,Q22)*pow(m_Mp*m_Mp/(2*nu(W2,-Q12,-Q22)),2)*
    (1+pow((2*Q12*Q22-nu(W2,-Q12,-Q22)*(Q12+Q22)),2)/(6*m_Mp*m_Mp*m_Mp*m_Mp*nu(W2,-Q12,-Q22)*nu(W2,-Q12,-Q22)));
}

double SigmaBW::fTT_2M(double W2,double Q12,double Q22)
{
  return kappa(W2,Q12,Q22)*pow(X(W2,Q12,Q22)/(nu(W2,-Q12,-Q22)*nu(W2,-Q12,-Q22)),3);
} 

// f_TS, J=0-, 0+, 1+, 2+, 2-
double SigmaBW::fTS_0M(double W2,double Q12,double Q22)
{
  return 0;
}

double SigmaBW::fTS_0P(double W2,double Q12,double Q22)
{
  return 0;
}

double SigmaBW::fTS_1P(double W2,double Q12,double Q22)
{
  return 2*kappa(W2,Q12,Q22)*m_Mp*m_Mp/(2*nu(W2,-Q12,-Q22))*Q22/(2*nu(W2,-Q12,-Q22))*
    (nu(W2,-Q12,-Q22)+Q12)/nu(W2,-Q12,-Q22)*(nu(W2,-Q12,-Q22)+Q12)/nu(W2,-Q12,-Q22);
}
  
double SigmaBW::fTS_2P(double W2,double Q12,double Q22)
{
  return  kappa(W2,Q12,Q22)*m_Mp*m_Mp*Q22*(nu(W2,-Q12,-Q22)-Q12)*(nu(W2,-Q12,-Q22)-Q12)/(4*pow(nu(W2,-Q12,-Q22),4));
}

double SigmaBW::fTS_2M(double W2,double Q12,double Q22)
{
  return 0;
} 

// f_ST, J=0-, 0+, 1+, 2+, 2-
double SigmaBW::fST_0M(double W2,double Q12,double Q22)
{
  return 0;
}

double SigmaBW::fST_0P(double W2,double Q12,double Q22)
{
  return 0;
}

double SigmaBW::fST_1P(double W2,double Q12,double Q22)
{
  return fTS_1P(W2,Q22,Q12);
}
  
double SigmaBW::fST_2P(double W2,double Q12,double Q22)
{
  return fTS_2P(W2,Q22,Q12);
}

double SigmaBW::fST_2M(double W2,double Q12,double Q22)
{
  return 0;
} 

// f_SS, J=0-, 0+, 1+, 2+, 2-
  
double SigmaBW::fSS_0M(double W2,double Q12,double Q22)
{
  return 0;
}

double SigmaBW::fSS_0P(double W2,double Q12,double Q22)
{
  return 2*kappa(W2,Q12,Q22)*pow(m_Mp*m_Mp*sqrt(Q12*Q22)/(3*nu(W2,-Q12,-Q22)*nu(W2,-Q12,-Q22)),2);
}

double SigmaBW::fSS_1P(double W2,double Q12,double Q22)
{
  return 0;
}
  
double SigmaBW::fSS_2P(double W2,double Q12,double Q22)
{
  return kappa(W2,Q12,Q22)*m_Mp*m_Mp*m_Mp*m_Mp*Q12*Q12*Q22*Q22/(3*pow(nu(W2,-Q12,-Q22),4));
}

double SigmaBW::fSS_2M(double W2,double Q12,double Q22)
{
  return 0;
}
//------------------------------------

//-------------------fAB------------------
//double SigmaBW::fAB(int a, int b, double W2,double Q12,double Q22, double Mp, int JP)
double SigmaBW::fAB(int a, int b, double W2,double Q12,double Q22)
{
  //Mp=M for all mesons except for pi0,eta and eta1, for which we take rho-meson mass!
  
  
  //1:jp0P, 2:jp0M, 3: jp1P, 4: jp2P, 5: jp2M
  
  if(m_JP==2){
    //if(a==1 && b==1){	return fTT_0M(W2,Q12,Q22,pow(Mp,2)); }
    if(a==1 && b==1){	return fTT_0M(W2,Q12,Q22); }
    else if(a==1 && b==0){ return fTS_0M(W2,Q12,Q22); }
    else if(a==0 && b==1){ return fST_0M(W2,Q12,Q22); }
    else if(a==0 && b==0){ return fSS_0M(W2,Q12,Q22); }
    else{  std::cout << "Error: a,b = 0 or 1" << std::endl; return 0;}
  }
  else if(m_JP==1){
    if(a==1 && b==1){	return fTT_0P(W2,Q12,Q22); }
    else if(a==1 && b==0){ return fTS_0P(W2,Q12,Q22); }
    else if(a==0 && b==1){ return fST_0P(W2,Q12,Q22); }
    else if(a==0 && b==0){ return fSS_0P(W2,Q12,Q22); }
    else{  std::cout << "Error: a,b = 0 or 1" << std::endl; return 0;}
  }
  else if(m_JP==3){
    if(a==1 && b==1){	return fTT_1P(W2,Q12,Q22); }
    else if(a==1 && b==0){ return fTS_1P(W2,Q12,Q22); }
    else if(a==0 && b==1){ return fST_1P(W2,Q12,Q22); }
    else if(a==0 && b==0){ return fSS_1P(W2,Q12,Q22); }
    else{  std::cout << "Error: a,b = 0 or 1" << std::endl; return 0;}
  }
  else if(m_JP==4){
    if(a==1 && b==1){	return fTT_2P(W2,Q12,Q22); }
    else if(a==1 && b==0){ return fTS_2P(W2,Q12,Q22); }
    else if(a==0 && b==1){ return fST_2P(W2,Q12,Q22); }
    else if(a==0 && b==0){ return fSS_2P(W2,Q12,Q22); }
    else{  std::cout << "Error: a,b = 0 or 1" << std::endl; return 0;}
  }
  else if(m_JP==5){
    if(a==1 && b==1){	return fTT_2M(W2,Q12,Q22); }
    else if(a==1 && b==0){ return fTS_2M(W2,Q12,Q22); }
    else if(a==0 && b==1){ return fST_2M(W2,Q12,Q22); }
    else if(a==0 && b==0){ return fSS_2M(W2,Q12,Q22); }
    else{  std::cout << "Error: a,b = 0 or 1" << std::endl; return 0;}
  }
  else{  std::cout << "SigmaBW.cxx, after fAB: This JP doesn't used" << std::endl; return 0;}
    
}
//-------------------------------------

// BW ------------------------
double SigmaBW::BW(double W2)
{
  return 1/M_PI*m_M*m_Gamma/((m_M*m_M-W2)*(m_M*m_M-W2)+m_M*m_M*m_Gamma*m_Gamma);
}
//-----------------------------

double SigmaBW::sigmaAB_BW(int a, int b,double W2,double Q12,double Q22)
{
  return (2*m_J+1)*8*M_PI*M_PI*GammaGG()/m_M*fAB(a,b,W2,Q12,Q22)*BW(W2);
}
//-----------------------------------------------

double SigmaBW::Sigma(double s,double s1,double s2,double t1,double t2,double W2,double phi,double theta)
{
  /*
  double M,Mp;
  int J,JP;
  double Gamma;
  */
/*  
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

  return 2*rho1PP*2*rho2PP*sigmaAB_BW(1,1,W2,-t1,-t2)+
    2*rho1PP*rho200*sigmaAB_BW(1,0,W2,-t1,-t2)+
    rho100*2*rho2PP*sigmaAB_BW(0,1,W2,-t1,-t2)+
    rho100*rho200*sigmaAB_BW(0,0,W2,-t1,-t2);
  //+ 2*Absrho1PM(s1,t1,t2,W2,m2)*Absrho2PM(s2,t1,t2,W2,m2)*tauTT(-t1,-t2,W2,theta)*cos2phi1(phi,s,s1,s2,W2,t1,t2)-
  //8*Absrho1P0(s1,t1,t2,W2,m2)*Absrho2P0(s2,t1,t2,W2,m2)*tauTS(-t1,-t2,W2,theta)*cosphi1(phi,s,s1,s2,W2,t1,t2);
  */
/*
    double  u1 = s2-me2-t1;
   double  u2 = s1-me2-t2;
   double  zeta_s2 = W2-t1-t2; //2*nu
   double  la = (t2-t1-W2)*(t2-t1-W2)/4. - t1*W2; // 1./4*lambda(t2,t1,W2)

   double rho1zz = (u2-zeta_s2/2)*(u2-zeta_s2/2)/la - 1;
   double rho2zz = (u1-zeta_s2/2)*(u1-zeta_s2/2)/la - 1;

   double rho1pp=rho1zz+2.+4.*me2/t1;
   double rho2pp=rho2zz+2.+4.*me2/t2;
   long double zero = 0.;
   double rho1pm = max(zero,((u2-zeta_s2/2)*(u2-zeta_s2/2)/la-1+4.*me2/t1)/2);
   double rho2pm = max(zero,((u1-zeta_s2/2)*(u1-zeta_s2/2)/la-1+4.*me2/t2)/2);
   double rho1pz = (u2-zeta_s2/2)*sqrt(rho1pm/la);
   double rho2pz = (u1-zeta_s2/2)*sqrt(rho2pm/la);

   long double sigma = ( rho1pp*rho2pp*sigmaTT(W2,t1,t2)
                         + rho1pp*rho2zz*sigmaTS(W2,t1,t2)
                         + rho1zz*rho2pp*sigmaST(W2,t1,t2)
                         + rho1zz*rho2zz*sigmaSS(W2,t1,t2)
                         //      - 8*rho1pz*rho2pz*tauTS*cos(phi)
                         //      + 2*rho1pm*rho2pm*tauTT*cos(2*phi)
                         );
*/
  return 1; //sigma;
}

