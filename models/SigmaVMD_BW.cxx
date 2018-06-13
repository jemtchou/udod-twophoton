
#define _USE_MATH_DEFINES //for M_PI

#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <numeric>
#include <cstdlib> // for exit(1)

#include "SigmaVMD_BW.hh"
#include "CkTable.h"
#include "Udod.h"
#include "UdodConstants.h"

using namespace std;
using namespace UDOD;

SigmaVMD_BW::SigmaVMD_BW(double M,int J,int P,double Gamma) 
{
  idchan =3;
  m_M = M;
  m_J = J;
  m_JP = P;
  m_Gamma = Gamma;
  mRhoSq=0.775*0.775;
};

double  SigmaVMD_BW::GammaGG()
{
  return 0.46;  

   //1:jp0P, 2:jp0M, 3: jp1P, 4: jp2P, 5: jp2M 
  if(m_JP==1){return 4;} 
  else if(m_JP==2){return 144;}
  else if(m_JP==3){return 32;}
  else if(m_JP==4){return 192./5;}
  else if(m_JP==5){return 64;}
  else{ cout << "SigmaVMD_BW.cxx: Error JP " << endl; return 0;}
}

// BW ------------------------
double  SigmaVMD_BW::BW(double W2)
{
  return 1;
  //return 1/M_PI*m_M*m_Gamma/((m_M*m_M-W2)*(m_M*m_M-W2)+m_M*m_M*m_Gamma*m_Gamma);
}


//a,b,Mp - only BW:
//double SigmaVMD_BW::sigmaAB_BW_or_VMD_BWf(int a, int b,double W2,double Q12,double Q22, double M, double Mp, int J, int JP,double Gamma) 
double SigmaVMD_BW::sigmaAB_VMD_BW(double W2,double Q12,double Q22) 
{
  return (2*m_J+1)*8*M_PI*M_PI*GammaGG()/m_M*(mRhoSq/(mRhoSq+Q12))*(mRhoSq/(mRhoSq+Q12))*(mRhoSq/(mRhoSq+Q22))*(mRhoSq/(mRhoSq+Q22))*BW(W2);
}
//-----------------------------------------------


double SigmaVMD_BW::Sigma(double s,double s1,double s2,double t1,double t2,double W2,double phi,double theta)
{
  /*
  double M;
  int J,JP;
  double Gamma;
  */
  
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
  
  double sigma = 2*rho1PP*2*rho2PP*sigmaAB_VMD_BW(W2,-t1,-t2)+
    2*rho1PP*rho200*sigmaAB_VMD_BW(W2,-t1,-t2)+
    rho100*2*rho2PP*sigmaAB_VMD_BW(W2,-t1,-t2)+
    rho100*rho200*sigmaAB_VMD_BW(W2,-t1,-t2);
  return sigma;
  cout << "sigma VMD BW: " << sigma << endl;
}
