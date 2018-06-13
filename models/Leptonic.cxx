#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <numeric>

#include "Leptonic.hh"
#include "UdodConstants.h"

using namespace std;
using namespace UDOD;


double Leptonic::q1muq2mu(double W2,double t1,double t2){
  return 1./2*(W2-t1-t2);
}

double Leptonic::x(double W2,double t1,double t2){
  return (q1muq2mu(W2,t1,t2)*q1muq2mu(W2,t1,t2)-t1*t2)/W2;
}

double Leptonic::Delta_t(double W2,double t1,double t2){
  return sqrt(4*x(W2,t1,t2)*(W2-4*m_m2));
}
double Leptonic::T(double W2,double t1,double t2){
  return 4*x(W2,t1,t2)*m_m2+t1*t2;
}
double Leptonic::L(double t1,double t2,double W2){
  return log(pow((q1muq2mu(W2,t1,t2)+1./2*Delta_t(W2,t1,t2)),2)/(4*x(W2,t1,t2)*m_m2+t1*t2));
}


//-------------sigmaAB--------------
double Leptonic::sigmaTT(double W2,double t1,double t2)
{
  long double q12 = 0.5*W2-0.5*t1-0.5*t2;
  long double x  = (q12*q12-t1*t2)/W2;
  long double T  = 4.*m_m2*x+t1*t2;
  long double deltat  = 2.*sqrt(x*(W2-4.*m_m2));
  long double L  = log((q12+0.5*deltat)*(q12+0.5*deltat)/(4.*m_m2*x+t1*t2));

  long double sig = pi*alpha*alpha/W2/x*(q12*L*(2.+2.*m_m2/x-4.*m_m2*m_m2/q12/q12+(t1+t2)/x+0.5*t1*t2*W2/x/q12/q12+0.75*t1*t1*t2*t2/x/x/q12/q12) - deltat*(1+m_m2/x+(t1+t2)/x +t1*t2/T+0.75*t1*t2/x/x));

  return sig; // GeV^-2
}

double Leptonic::sigmaTS(double W2,double t1,double t2){
  double q12 = 0.5*W2-0.5*t1-0.5*t2;
  double x  = (q12*q12-t1*t2)/W2;
  double T  = 4.*m_m2*x+t1*t2;
  double deltat  = 2.*sqrt(x*(W2-4.*m_m2));
  double L  = log((q12+0.5*deltat)*(q12+0.5*deltat)/(4.*m_m2*x+t1*t2));

  double sig = -pi*alpha*alpha*t2/W2/x/x*(deltat*(1.+t1/T*(6.*m_m2+t1+1.5*t1*t2/x))-L/q12*(4.*m_m2*x+t1*(W2+2.*m_m2)+t1*(t1+t2+1.5*t1*t2/x)));

  return sig;
} 


double Leptonic::sigmaST(double W2,double t1,double t2){
  return  sigmaTS(W2,t2,t1);
}

double Leptonic::sigmaSS(double W2,double t1,double t2){
  double q12 = 0.5*W2-0.5*t1-0.5*t2;
  double x  = (q12*q12-t1*t2)/W2;
  double T  = 4.*m_m2*x+t1*t2;
  double deltat  = 2.*sqrt(x*(W2-4.*m_m2));
  double L  = log((q12+0.5*deltat)*(q12+0.5*deltat)/(4.*m_m2*x+t1*t2));

  double sig = pi*alpha*alpha*t1*t2/W2/x/x/x*(L/q12*(2.*W2*x+3.*t1*t2)-deltat*(2.+t1*t2/T));

  return sig;
}

double Leptonic::tauTT(double W2,double t1,double t2){
  double q12 = 0.5*W2-0.5*t1-0.5*t2;
  double x  = (q12*q12-t1*t2)/W2;
  double deltat  = 2.*sqrt(x*(W2-4.*m_m2));
  double L  = log((q12+0.5*deltat)*(q12+0.5*deltat)/(4.*m_m2*x+t1*t2));

  double sig = -0.25*pi*alpha*alpha/W2/x*(2.*deltat/x*(2.*m_m2+(t1-t2)*(t1-t2)/W2+1.5*t1*t2/x)+L/q12*(16.*m_m2*m_m2-16.*m_m2*(t1+t2)-4.*t1*t2*(2.+2.*m_m2/x+(t1+t2)/x+0.75*t1*t2/x/x)));

  return sig;
}

double Leptonic::tauTS(double W2,double t1,double t2){
  double q12 = 0.5*W2-0.5*t1-0.5*t2;
  double x  = (q12*q12-t1*t2)/W2;
  double deltat  = 2.*sqrt(x*(W2-4.*m_m2));
  double L  = log((q12+0.5*deltat)*(q12+0.5*deltat)/(4.*m_m2*x+t1*t2));

  double sig = -pi*alpha*alpha*sqrt(t1*t2)/W2/x/x*(L*(2.*m_m2+t1+t2+1.5*t1*t2/x)+deltat*(2.-1.5*q12/x));

  return sig;
}

double Leptonic::Sigma(double s,double s1,double s2,double t1,double t2,double W2,double cphi,double theta)
{
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

   double cosphi = cphi;
   double cos2phi = 2*cphi*cphi-1;

   long double sigma = ( rho1pp*rho2pp*sigmaTT(W2,t1,t2)
                         + rho1pp*rho2zz*sigmaTS(W2,t1,t2)
                         + rho1zz*rho2pp*sigmaST(W2,t1,t2)
                         + rho1zz*rho2zz*sigmaSS(W2,t1,t2)
                         - 8*rho1pz*rho2pz*tauTS(W2,t1,t2)*cosphi
                         + 2*rho1pm*rho2pm*tauTT(W2,t1,t2)*cos2phi
                         );

  return sigma;
}

// Budnev, p.276
double Leptonic::SigmaTheorEE(double s,double s1,double s2,double t1,double t2,double W2,double phi,double theta)
{
  double p1sq=me2-(sqrt(s1)-t2)*(sqrt(s1)-t2);
  double p2sq=me2-(sqrt(s2)-t1)*(sqrt(s2)-t1);
  
  double p1mup2mu = s-p1sq-p2sq;
  double l= log(2*p1mup2mu/me2);

  
  double A = 6.36;
  double B = 15.7;
  double C = -13.8;

  return 28*pow(alpha,4)/(27*M_PI*me2)*(l*l*l-A*l*l+B*l+C)*
  nbarn;
  //15.396*pow(10,6)/(4*M_PI*M_PI);//nBn;
}


// Budnev, p.277
double Leptonic::SigmaTheorMuMu(double s,double s1,double s2,double t1,double t2,double W2,double phi,double theta)
{
  //double lmu = log(mmu2/me2); //10.7;
  double lmu = 10.7;
  double p1sq=me2-(sqrt(s1)-t2)*(sqrt(s1)-t2);
  double p2sq=me2-(sqrt(s2)-t1)*(sqrt(s2)-t1);
  
  double p1mup2mu = s-p1sq-p2sq;
  double l= log(2*p1mup2mu/me2);
  /*
   cout << "p1sq: " << p1sq << endl;
   cout << "p2sq: " << p2sq << endl;
   cout << "p1mup2mu: " << p1mup2mu << endl;
   cout << "l: " << l << endl;
  */
  
  double A = 6.36;
  double B = -258.;
  double C = 1850.;

  return 28*pow(alpha,4)/(27*M_PI*mmu2)*(l*l*l-A*l*l+B*l+C)*
    nbarn;
  //15.396*pow(10,6)/(4*M_PI*M_PI); //nBn
}

//Budnev, p.184

double Leptonic::SigmaTheorll(double s,double s1,double s2,double t1,double t2,double W2,double phi,double theta)
{
  return 28*pow(alpha,4)/(27*M_PI*m_m2)*log(s/me2)*log(s/me2)*log(s/m_m2)*
  nbarn;
  //15.396*pow(10,6)/(4*M_PI*M_PI); //nBn
}

//Budnev, p.193
double Leptonic::SigmaTheorEE1(double s,double s1,double s2,double t1,double t2,double W2,double phi,double theta)
{
  double p1sq=me2-(sqrt(s1)-t2)*(sqrt(s1)-t2);
  double p2sq=me2-(sqrt(s2)-t1)*(sqrt(s2)-t1);
  
  double p1mup2mu = s-p1sq-p2sq;
  double l= log(2*p1mup2mu/me2);
  
  return 28*pow(alpha,4)/(27*M_PI*me2)*l*l*l*
    nbarn;//nbn
  //15.396*pow(10,6)/(4*M_PI*M_PI); //nBn
}
