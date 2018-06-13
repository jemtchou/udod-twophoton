// version 2.0 based on Galuga
#include <cmath>
#include <cfloat>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <iomanip>

//#include "CLHEP/Random/Ranlux64Engine.h"
//#include "CLHEP/Random/RandFlat.h"

#include "cuba_int.h"

#include "Udod.h"
//#include "UdodChannel.h"
#include "UdodChannelRegistry.h"

#include "UdodConstants.h"
#include "Udod.icc"

using namespace std;
using namespace UDOD;

//------------------------------------------------------------------------------
Udod::Udod(double Ecm)
{
   // set beam energy
   rs = Ecm;
   s = rs*rs;
   beta2= 1.-4.*me2/s;
   beta = sqrt(beta2);

   // set default cuts
   fEnergyCut[0].SetCut(me,rs/2,"e- energy","GeV");
   fEnergyCut[1].SetCut(me,rs/2,"e+ energy","GeV");
   fThetaCut[0].SetCut(0.,180.,"e- theta","deg");
   fThetaCut[1].SetCut(0.,180.,"e+ theta","deg");
   fXMassCut.SetCut(m_pi, rs-2.*me,"Hadron mass","GeV");

   // set cross-section extremum
   fSigmaMax=0;
   fSigmaMin=DBL_MAX;

   fGammaModel = GVMD;

   fChannelRegistry = new UdodChannelRegistry(this);

   fEventsTotal = 0;
   fEventsAccepted = 0;
   fVerbosity = 3;

   fChannel = 0;
   DFS.reserve(24);
}

//------------------------------------------------------------------------------
Udod::~Udod()
{
}

void Udod::RegisterChannel(UdodChannelBase* ch)
{
  fChannelRegistry->AddChannel(ch);
}

//------------------------------------------------------------------------------
void Udod::Message(string mess, unsigned int level)
{
   if( level > 4 ) {
      level = 4;
   }
   if( level >= fVerbosity ) {
      cout << DebugLevel[level] << ": " << mess << endl;
   }
   cout.flush();
}

//------------------------------------------------------------------------------
// random number from zero to one
double Udod::Random()
{
  return fEngine.Random();
}

//------------------------------------------------------------------------------
// array of random numbers form 0 to 1
void Udod::RandomArray(int size, double* rndarray)
{
  fEngine.RandomArray(size, rndarray);
}

//------------------------------------------------------------------------------
void Udod::SetEnergyCut(int flag, double min, double max)
{
   double mintmp = min;
   double maxtmp = max;

   if( mintmp > maxtmp ) {
      Message("Wrong cut on electron/positron energy: max < min."
              " Ignore this cut",ERROR);
      return;
   }

   if( mintmp < me ) {
      mintmp = me;
   }
   if( maxtmp > rs/2 ) {
      maxtmp = rs/2;
   }

   fEnergyCut[flag].SetCut(mintmp,maxtmp,ParticleName[flag]+" energy","GeV");
}

//------------------------------------------------------------------------------
void Udod::SetThetaCut(int flag, double min, double max)
{
   double mintmp = min;
   double maxtmp = max;

   if( mintmp > maxtmp ) {
      Message("Wrong polar angle cut: max < min. Ignore this cut",ERROR);
      return;
   }

   if( mintmp < 0.) {
      mintmp = 0.;
   }
   if( maxtmp > 180.) {
      maxtmp = 180.;
   }

   fThetaCut[flag].SetCut(mintmp*pi/180.,maxtmp*pi/180,ParticleName[flag]+" theta","rad");
}

//------------------------------------------------------------------------------
void Udod::SetHadronMassCut(double min, double max)
{
   double mintmp = min;
   double maxtmp = max;
   if( mintmp > maxtmp ) {
      Message("Wrong hadron mass cut: max < min. Ignore this cut",ERROR);
      return;
   }

   if( mintmp < m_pi ) {
      mintmp = m_pi;   // minimal hadron mass - pion mass
   }
   if( maxtmp > rs-2.*me ) {
     maxtmp = rs-2.*me;   // maximum of hadron mass
   }
   
   fXMassCut.SetCut(mintmp,maxtmp,"Hadron mass","GeV");
}

//------------------------------------------------------------------------------
bool Udod::Initialize()
{
   fThetaCut[0].Print();
   fThetaCut[1].Print();
   fEnergyCut[0].Print();
   fEnergyCut[1].Print();
   fXMassCut.Print();

   cout << endl << " ------- CUBA ------- " << endl;
   cuba_int cuba;
   cuba.set_fun( fastdelegate::MakeDelegate(this, &Udod::Integrand) );

   // set parameters the same as in 'vegas()' function
   int ncomp = fChannelRegistry->GetNChannels();
   cuba.set_dim(5);
   cuba.set_ncomp(ncomp);
   cuba.set_nstart(2000000);
   cuba.set_nincrease(1000);
   cuba.set_nom(VEGAS);
   cuba.set_verbose(1);

   cuba.Integral();

   cout << setprecision(8) << " " << cuba.str_nom() <<" RESULT " << endl;
   int i;
   for (i=0; i<ncomp; i++)
     {
       fChannelRegistry->GetChannel(i)->SetIntegralCS(cuba.get_integral(i));

       cout << "Channel " << fChannelRegistry->GetChannel(i)->GetName() << endl;
       cout << cuba.get_integral(i) << " [nb] +/- " << cuba.get_error(i)
	    << " Prob " << setprecision(3) << cuba.get_prob(i) << endl;
       cout <<  setprecision(18) << " Max. cross-section: " << fChannelRegistry->GetChannel(i)->GetMaximumCS() << setprecision(3) << endl;
     }
   
   fChannelRegistry->Initialize();

   return true;
}

//------------------------------------------------------------------------------
void Udod::Finalize()
{
   Message("Finalize");
   cout << "Total attempts: " << fEventsTotal
        << " accepted " << fEventsAccepted << endl;
}

//------------------------------------------------------------------------------
long double Udod::ComputePhaseSpace()
{
   double rand[5];
   RandomArray(5, rand);

   //------------ TEST
   // rand[0]=2.145740389823916E-002;
   // rand[1]=4.291901985804243E-002;
   // rand[2]=2.985138694445294E-002;
   // rand[3]=3.704054156939191E-002;
   // rand[4]=0.240349367260933;
   //-----------------
  
   phsp = 0;
   W2 = 0.;
   s1 = 0.; s2 = 0.; t1 = 0.; t2 = 0.;
   ds1 = 0.; ds2 = 0.; dt1 = 0.; dt2 = 0.;

   // calculate mass of hadron state
   double WMin = fXMassCut.Min();
   double WMax = fXMassCut.Max();
     
   // logarithmic mapping
   dtau   = 2*logl(WMax/WMin);
   W2= WMin*WMin*exp(dtau*rand[4]);
   W = sqrt(W2);
  
   if(! CalculateT2(rand[0]))
     return 0;
   //tCut[1].Print();
   // cout << "t2 = " << t2 << " " << dt2 << endl;
   //UdodHisto::GetInstance()->Fill(4,-t2,1);
  
   if(! CalculateT1(rand[1])) // using t2 as a kinematic limit
      return 0;
   //tCut[0].Print();
   //cout << "t1 = " << t1 << " " << dt1 << endl;
   //UdodHisto::GetInstance()->Fill(3,-t1,1);

   if(! CalculateS1(rand[2])) // using limits from t1, t2
     return 0;
   //sCut[0].Print();
   //cout << "s1 = " << s1 << " " << ds1 << endl; 
   //UdodHisto::GetInstance()->Fill(1,s1,1);


   if(! CalculateS2(rand[3])) // using limits from t1, t2, s1
     return 0;
   //sCut[1].Print();
   //cout << "s2 = " << s2 << " " << ds2 << endl;
   //UdodHisto::GetInstance()->Fill(2,s2,1);

   phsp = dt1*dt2*ds1*ds2*dtau*W2/s*sqrtl(lambda(W2,t1,t2))*alpha*alpha/(16*pi*pi*beta2*s);

   CalculateMomenta();

   // phsp =  dt1*dt2*ds1*ds2*dtau*W2/s*pi*pi*t1*t2/beta2/s;
        long double test[4];
   memset(test,0,sizeof(test));
   int ipart = 0;
   for(ipart = 0; ipart < 2; ipart++) {
      ComputeElectronScattering(ipart);
   }
   ComputeGammaGammaState();

   return phsp;
}

bool Udod::CalculateT2(const double rand)
{
  long double s1 = 0.;
  // t2max
  long double s1u = (rs-me)*(rs-me);
  long double s1l = (W+me)*(W+me);
  s1u = std::min(s1u, me2 + s*(1 - 2*fEnergyCut[1].Min()/rs));
  s1l = std::max(s1l, me2 + s*(1 - 2*fEnergyCut[1].Max()/rs));
  if (s1l>=s1u) 
    return false;
  
  double theta =  fThetaCut[1].Min();
  long double tmp = 4*me2/s + beta2*sinl(theta)*sinl(theta);
  long double s1hat = me2 + s*beta2*sinl(theta)*sinl(theta)/(tmp*(1 + 2*me/sqrtl(tmp*s)));
  long double t2hat = 2*me2 - me*rs*sqrtl(4*me2/s + beta2*sinl(theta)*sinl(theta));
  
  long double t2max = 0.;
  if(s1hat >= s1u)
    {
      s1 = s1u;
      if( fEnergyCut[1].Min() <= me )
	{
	   tmp = 0.;
	   t2max = 2*me2 - me*rs;
	}
      else
	{
	  tmp = lambda(s, s1, me2);
	  t2max = (3*me2 -s + s1 -cosl(theta)*beta*sqrtl(tmp))/2;
	}
      t2max = (me2*(s1 - me2)*(s1 - me2)/s + beta2*tmp*sinl(theta)*sinl(theta)/4.)/t2max;
    }
  else if(s1hat <= s1l )
     {
       s1 = s1l;
       tmp = lambda(s, s1, me2);
       t2max = (3*me2 -s + s1 -cosl(theta)*beta*sqrtl(tmp))/2;
       t2max = (me2*(s1 - me2)*(s1 - me2)/s + beta2*tmp*sinl(theta)*sinl(theta)/4.)/t2max;
     }
  else
    {
       t2max = t2hat;
    }
  
  // t2min
  theta =  fThetaCut[1].Max();
  s1 = s1l;
  tmp = lambda(s,s1,me2);
  long double t2min = (3*me2 - s + s1 + cosl(theta)*beta*sqrtl(tmp))/2;

  s1 = s1u;
  long double t2l = 0;
  if(fEnergyCut[1].Min() <= me )
    {
      tmp = 0.;
      t2l = 2*me2 - me*rs;
     }
  else
    {
      tmp = lambda(s, s1, me2);
      t2l = (3*me2 - s + s1 + cosl(theta)*beta*sqrtl(tmp))/2;
    }
  t2min = std::min(t2min, t2l);
  long double t2umin = -1e12;
  t2min = std::max(t2min, t2umin);
  tCut[1].SetCut(t2min,t2max,"t2","GeV^2");
  
  dt2 = logl(t2max/t2min);
  t2 = t2min*expl(dt2*rand);

  return true;
}

bool Udod::CalculateT1(const double rand)
{
  // t1max
  long double s2 = 0.;
  long double s2u = (rs - me)*(rs - me);
  long double s2l = (W + me)*(W + me);
  s2u = std::min(s2u, me2 + s*(1 - 2*fEnergyCut[0].Min()/rs));
  s2l = std::max(s2l, me2 + s*(1 - 2*fEnergyCut[0].Max()/rs));
  
  double theta =  fThetaCut[0].Min();
  long double tmp = 4*me2/s + beta2*sinl(theta)*sinl(theta);
  long double s2hat = me2 + s*beta2*sinl(theta)*sinl(theta)/(tmp*(1 + 2*me/sqrtl(tmp*s)));
  long double t1hat = 2*me2 - me*rs*sqrtl(4*me2/s + beta2*sinl(theta)*sinl(theta));
 
  long double t1max = 0.;
  if(s2hat >= s2u)
    {
      s2 = s2u;
      if( fEnergyCut[0].Min() <= me )
	{
	  tmp = 0.;
	  t1max = 2*me2 - me*rs;
	}
      else
	{
	  tmp = lambda(s, s2, me2);
	  t1max = (3*me2 -s + s2 - cosl(theta)*beta*sqrtl(tmp))/2;
	}
      t1max = (me2*(s2 - me2)*(s2 - me2)/s + beta2*tmp*sinl(theta)*sinl(theta)/4.)/t1max;
    }
  else if(s2hat <= s2l )
    {
      s2 = s2l;
      tmp = lambda(s, s2, me2);
      t1max = (3*me2 - s + s2 -cosl(theta)*beta*sqrtl(tmp))/2;
      t1max = (me2*(s2 - me2)*(s2 - me2)/s + beta2*tmp*sinl(theta)*sinl(theta)/4.)/t1max;
    }
  else
    {
      t1max = t1hat;
    }

  // t1min
  theta =  fThetaCut[0].Max();
  s2 = s2l;
  tmp = lambda(s,s2,me2);
  long double t1min = (3*me2 - s + s2 + cosl(theta)*beta*sqrtl(tmp))/2;

  s2 = s2u;
  long double t1l = 0;
  if(fEnergyCut[0].Min() <= me )
    {
      tmp = 0.;
      t1l = 2*me2 - me*rs;
     }
  else
    {
      tmp = lambda(s, s2, me2);
      t1l = (3*me2 - s + s2 + cosl(theta)*beta*sqrtl(tmp))/2;
    }
  t1min = std::min(t1min, t1l);
  
  // Kinematical limits due to t2
  long double y2 = sqrtl(1. - 4*me2/t2);
  long double Q = 4*(s + t2 - 4*me2)/(1 + beta*y2) - t2 - W2;
  long double a1 = 2*(Q + t2 + 2*me2 + W2);
  long double b1 = Q*Q -W2*W2 + 2*W2*t2 - t2*t2 -8*me2*t2 - 8*me2*W2;
  long double c1 = 4*me2*(W2 - t2)*(W2 - t2);
  long double D1 = (Q + t2 - W2 + 4*me*W)*(Q + t2 - W2 - 4*me*W)*(Q*Q - 2*Q*t2 + 2*Q*W2 + t2*t2 + W2*W2 - 16*me2*t2 - 2*W2*t2);  

  long double tmp1 = sqrtl(D1)/a1;
  long double t1low = -(b1/a1 + tmp1)/2;
  long double t1upp = c1/(a1*t1low);

  t1min = std::max(t1min, t1low);
  t1max = std::min(t1max, t1upp);
  if(t1min >= t1max)
    return false;
  tCut[0].SetCut(t1min,t1max,"t1","GeV^2");

  dt1 = logl(t1max/t1min);
  t1 = t1min*expl(dt1*rand);

  return true;
}

bool Udod::CalculateS1(const double rand)
{
  long double y1 = sqrtl(1-4*me2/t1);
  long double y2 = sqrtl(1-4*me2/t2);
  long double KW = 0.5*sqrtl(lambda(W2,t1,t2));
  long double nu = (W2 - t1 - t2)/2;
  
  long double s1min = me2 + (W2 - t1 + t2 + 2*y1*KW)/2;
  long double s1max = me2 + 2*(s + t2 -4*me2)/(1 + beta*y2);
  sCut[0].SetCut(s1min, s1max, "s1","GeV^2");

  ds1 = logl(s*(1 + beta)*(1+beta)/((nu + KW)*(1 + y1)*(1 + y2)));
  long double X = (nu + KW)*(1 + y1)*expl(ds1*rand);

  s1 = X/2 + me2 + t2 + 2*me2*t2/X; 
  UdodHisto::GetInstance()->Fill(20, s1min,1);
  UdodHisto::GetInstance()->Fill(21, s1max,1);

  return true;
}

bool Udod::CalculateS2(const double rand)
{
  // hep-ex/9710506, formulae at pp.10-11
  long double A = lambda(s1, t2, me2);
  long double B = -2*s*me2*t1 - 2*me2*s1*s1 + 8*t2*me2*me2 - 2*me2*t2*t2 - 2*s*s1*W2 + 2*me2*s*W2 + 2*t1*s*s1 + 2*s*t2*s1 + 4*me2*s1*W2 + 4*me2*me2*s1 + 2*t1*t2*s - 2*t2*me2*t1 - 2*t2*t2*s - 2*me2*t2*s + 2*t1*t2*s1 - 4*me2*me2*W2 - 2*t1*s1*s1 + 2*s*t2*W2 - 2*me2*me2*me2 + 2*me2*me2*t1;  
  long double C = -2*s*me2*me2*W2 - 2*t1*t1*me2*s1 - 2*t1*t2*s*s + 2*t1*t2*s*s1 - 2*t1*t1*s*s1 + t1*t1*s*s + t1*t1*s1*s1 + t2*t2*s*s + me2*me2*s1*s1 + me2*me2*t1*t1 - 6*me2*me2*t1*me2 - 2*me2*me2*s1*me2 - 4*me2*me2*s1*W2 + 2*me2*me2*t2*s + 2*me2*me2*t2*s1 + 8*me2*me2*t1*s1 - 2*t2*s*s*W2 - 2*t1*s*s*W2 - 2*me2*t1*s1*s1 +me2*me2*me2*me2 - 2*me2*s*t2*s1 + 4*me2*me2*me2*W2 + me2*me2*t2*t2 + 4*me2*t1*t2*s - 2*me2*t1*t2*s1 -6*me2*me2*me2*t2 +s*s*W2*W2 + 6*me2*s*t2*W2 - 4*s*me2*W2*W2 - 2*s*me2*t1*t1 + 2*s*t1*me2*me2 -2*s*me2*t2*t2 +2*t1*t2*me2*me2 + 2*s*me2*s1*W2 -4*s*t1*t2*W2 - 2*s*t1*me2*s1 +6*s*t1*me2*W2 + 2*s*t1*s1*W2;     

  long double G3 = me2*s1*s1 - 2*me2*me2*s1 - s*t2*s1 - 3*me2*t2*s + me2*me2*me2 + t2*s*s + t2*t2*s;
  long double G4 = -2*t1*me2*s1 - t2*t1*me2 +me2*me2*t1 - me2*W2*t1 + me2*t2*t2 + t1*t2*W2 - t1*s1*W2 - 2*me2*t2*W2 + me2*W2*W2 + t1*s1*s1 + t1*t1*s1 - t1*t2*s1;
  long double D =  16*G3*G4;
  if(D<0) return false;

  long double s2max = (sqrtl(D) - B)/2/A;
  long double s2min = C/(A*s2max);
  sCut[1].SetCut(s2min, s2max, "s2","GeV^2");

  s2 = (sqrtl(D)*sinl((rand-0.5)*pi) - B)/2/A;
  ds2 = 1;
  UdodHisto::GetInstance()->Fill(22, s2min,1);
  UdodHisto::GetInstance()->Fill(23, s2max,1);
  
  return true;
}

bool Udod::CalculateMomenta()
{
  E1 = (s + me2 - s2)/2/rs;
  E2 = (s + me2 - s1)/2/rs;
  EX = (s1 + s2 -2*me2)/2/rs;
  P1 = sqrtl(lambda(s,s2,me2));
  P2 = sqrtl(lambda(s,s1,me2));
  PX = sqrtl(EX*EX - W2);

  long double G3 = me2*s1*s1 - 2*me2*me2*s1 - s*t2*s1 - 3*me2*t2*s + me2*me2*me2 + t2*s*s + t2*t2*s;
  long double G4 = -2*t1*me2*s1 - t2*t1*me2 + me2*me2*t1 - me2*W2*t1 + me2*t2*t2 + t1*t2*W2 - t1*s1*W2 - 2*me2*t2*W2 + me2*W2*W2 + t1*s1*s1 + t1*t1*s1 - t1*t2*s1;
  long double G5 = -2*t2*me2*s2 - t1*t2*me2 + me2*me2*t2 - me2*W2*t2 + me2*t1*t1 + t2*t1*W2 - t2*s2*W2 - 2*me2*t1*W2 + me2*W2*W2 + t2*s2*s2 + t2*t2*s2 - t2*t1*s2;
  long double D1 = -(me2*s2*s2 - 2*me2*me2*s2 - s*t1*s2 - 3*me2*t1*s + me2*me2*me2 + t1*s*s + t1*t1*s)/4;
  long double D2 = -G5/4;
  long double D3 = -G3/4;
  long double D4 = -G4/4;
  long double D6 = s/8*(-(s - 4*me2)*(W2 - t1 -t2) + (s1 - t2 -me2)*(s2 - t1 -me2) + t1*t2) - me2/4*(s1 - me2)*(s2 - me2);
  long double D5 = D1 + D3 + 2*D6;
  long double D7 = 1./16*(2*W2*(s1*s2 - s*W2) - 2*t1*(-t1*s1 + s*t1 +s1*s2 +W2*s1 - 2*s*W2) - 2*t2*(-t2*s2 + s*t2 +s1*s2 +W2*s2 - 2*s*W2) +2*t1*t2*(-s1 - s2 + 2*s + 2*W2) - 2*me2*(me2*t2 - t2*t2 - t1*t1 - me2*W2 + me2*t1 - 2*W2*W2 + W2*s1 + 2*t1*t2 + 3*t1*W2 - t1*s1 + 3*t2*W2 - t2*s2 - t2*s1 + s2*W2 - t1*s2));
  
  sth1 = 2*sqrtl(D1)/s/beta/P1;
  sth2 = 2*sqrtl(D3)/s/beta/P2;
  sthX = 2*sqrtl(D5)/s/beta/PX;

  cospht = D7/(sqrtl(D2*D4));
  sinpht = sin(acos(cospht));

  return true;
}

//------------------------------------------------------------------------------
bool Udod::ComputeLimits(double W, UdodCut& tCut, UdodCut& sCut,
                         UdodCut& thetaCut, UdodCut& energyCut)
{
   double smin = max((me+W)*(me+W), me2 + s*(1.-2.*energyCut.Max()/rs));
   double smax = min( (rs-me)*(rs-me), me2 + s*(1.-2.*energyCut.Min()/rs));
   if( smin > smax) {
      return false;
   }

   long double sqlasi0=sqrtl(lambda(s,me2,smin));
   long double sqlasi1=sqrtl(lambda(s,me2,smax));

   double tmin = 1./2.*(-s+smin+3*me2+beta*sqlasi0*cos(thetaCut.Max()*to_rad) );
   long double sintheta = sin(thetaCut.Min()*to_rad);
   long double sinhalftheta = sin(thetaCut.Min()*to_rad/2);
   long double xi = 4*me2/s+beta2*sintheta;
   long double ZeroDerSPoint = me2 + s*beta2*sintheta*sintheta/xi/
                               (1+2*me/sqrt(s*xi));

   double tmax = 0;
   if ( ZeroDerSPoint <= smin ) {
      tmax = -2*me2*(me2-smin)*(me2-smin)/s/(beta*sqlasi0+s-smin-3*me2) -
             beta*sqlasi0*sinhalftheta*sinhalftheta;
   }

   if ( ZeroDerSPoint >= smax ) {
      tmax = -2*me2*(me2-smax)*(me2-smax)/s/(beta*sqlasi1+s-smax-3*me2) -
             beta*sqlasi1*sinhalftheta*sinhalftheta;
   }

   if ( ZeroDerSPoint< smax && ZeroDerSPoint> smin ) {
      tmax = 2*me2-me*rs*sqrt(xi);
   }

   if(tmin > tmax) {
      Message("tmax < tmin", FATAL);
      exit(1);
   }

   sCut.SetCut(smin,smax);
   tCut.SetCut(tmin,tmax);

   return true;
}

//------------------------------------------------------------------------------
void Udod::ComputeElectronScattering(int particle)
{
   long double s_1, t_1, s_2, t_2;
   int sign = 1;
   // int addpi = 0;
   if(particle == ELECTRON) {
      s_1 = s1; s_2 = s2; t_1 = t1; t_2 = t2;
   } else if (particle == POSITRON) {
      s_1 = s2; s_2 = s1; t_1 = t2; t_2 = t1;
      sign = -1;
   } else {
      Message("Requested particle neither e+ nor e-",FATAL);
      return;
   }

   long double E = (s - s_2 + me2)/2/rs;
   long double P = sqrtl(lambda(s,me2,s_2))/2/rs;
//   long double sint = ;
//  long double cost = sqrtl((1-sint)(1+sint))
//   long double Theta = acosl(sign*(s-s_2+2*t_1-3*me2)/(2*beta*rs*P));
//   long double Theta = acosl(cost);

   long double  Theta  = 2.*asinl(sqrtl(-t_1/(s-s_2)));
//	               th2    = 2d0*dAsin(dsqrt(-t2/(s-s2)))

   long double cosphi = coefGphi1(W2, s, s_1, t_1, s_2, t_2)/
                        sqrtl(coefG1(s, s_1, t_1, s_2, t_2)*
                              coefGX(W2, s, s_1, t_1, s_2, t_2));
   long double Phi = sign*acosl(cosphi);

   FS[particle].Set(E,P,Theta,Phi);
}

//------------------------------------------------------------------------------
void Udod::ComputeGammaGammaState()
{
   long double E = (s1+s2-2*me2)/2/rs;
   long double P = sqrtl((E-W)*(E+W));
   long double Theta = acosl((s2-s1+2*(t2-t1))/2/beta/rs/P);
   long double Phi = 0.;
   FS[XSTATE].Set(E,P,Theta,Phi);
}

//------------------------------------------------------------------------------
long double Udod::ComputeCrossSection()
{
   long double sigma;
   sigma = fChannel->Sigma(s,s1,s2,t1,t2,W2,cospht,sthX);
   sigma = sigma*phsp;

   UdodHisto::GetInstance()->Fill(5,sigma*nbarn,1);
   return sigma*nbarn; // GeV->nbarn
}

//------------------------------------------------------------------------------
double Udod::sigma_ab(int a, int b)
{
   // return sigma_s or sigma_t
   double sigma_ab = 0;

   if( fGammaModel == GVMD ) {
      sigma_ab = ht_or_s_GVMD(-t1,a)*ht_or_s_GVMD(-t2,b);
   } else if( fGammaModel == VMD ) {
      sigma_ab = ht_or_s_VMD (-t1,a)*ht_or_s_VMD (-t2,b);
   } else if( fGammaModel == RHOPOLE ) {
      sigma_ab = ht_or_s_pole(-t1,a)*ht_or_s_pole(-t2,b);
   } else {
      Message("No such GammaModel implemented", FATAL);
      exit(1);
   }

   return sigma_ab*1;
}

//------------------------------------------------------------------------------
double Udod::ht_or_s_GVMD(double Q2,int transversity)
{
   double r=3./4., ksi=1./4., m12=0.54, m22=1.8;
   double pQ_1=1+Q2/m12,pQ_2=1+Q2/m22;

   if( transversity==1 ) {
      return (r/pQ_1/pQ_1+(1.-r)/pQ_2);
   }
   if( transversity==0 ) {
      return (ksi*( r*Q2/m12/pQ_1/pQ_1+(1.-r)*(m22/Q2*log(pQ_2)-1./pQ_2) ));
   }
   cout << "Error in ht_or_s" << "   transversity=" << transversity << endl;
   exit(1);
}

//------------------------------------------------------------------------------
double Udod::ht_or_s_VMD(double Q2,int transversity)
{
   //!!!!!!!!!!! proverit' konstanti
   double r_rho=0.65, r_omega=0.08, r_phi=0.05, r_c=1-r_rho-r_omega-r_phi;
   double m_rho2=0.5920, m_omega2=0.611, m_phi2=1.039, m_02=1.96;

   if( transversity==1 ) {
      return ( r_rho/(1+Q2/m_rho2)/(1+Q2/m_rho2) +
               r_omega/(1+Q2/m_omega2)/(1+Q2/m_omega2) +
               r_phi/(1+Q2/m_phi2)/(1+Q2/m_phi2) +
               r_c/(1+Q2/m_02)
             );
   }

   if( transversity==0 ) {
      return ( r_rho/(1+Q2/m_rho2)/(1+Q2/m_rho2)*Q2/m_rho2/4 +
               r_omega/(1+Q2/m_omega2)/(1+Q2/m_omega2)*Q2/m_omega2/4 +
               r_phi/(1+Q2/m_phi2)/(1+Q2/m_phi2)*Q2/m_phi2/4
             );
   }
   cout << "Error in ht_or_s" << "   transversity=" << transversity << endl;
   exit(1);
}

//------------------------------------------------------------------------------
double Udod::ht_or_s_pole(double Q2,int transversity)
{
   //!!!!!!!!!!! proverit' konstanti
   double m_rho2=0.5920;

   double tmp = 1./(1+Q2/m_rho2)/(1+Q2/m_rho2);

   if( transversity==1 ) {
      return tmp;
   }
   if( transversity==0 ) {
      return tmp*Q2/m_rho2/4;
   }
   Message("Wrong transversity",FATAL);
   return 0;
}

//------------------------------------------------------------------------------
void Udod::Integrand(const double xx[], double* ff)
//------------------------------------------------------------------------------
{
   long double phsp = ComputePhaseSpace();
   int nchan = fChannelRegistry->GetNChannels();
   int ichan;
   for(ichan = 0; ichan<nchan; ichan++)
     {
       fChannel = fChannelRegistry->GetChannel(ichan);
       double sigma = ComputeCrossSection();

       if( sigma < 0 ) {
	 cerr << "ERROR: Negative TPA cross-section!" << sigma << endl;
	 cout << " t1=" << t1 << endl;
	 tCut[0].Print();
	 cout << " t2=" << t2 << endl;
	 tCut[1].Print();
	 cout << " W=" << W << endl;
	 fXMassCut.Print();
	 cout.flush();
	 exit(1);
       }
       
       // max and min values of cross section
       if( sigma > fChannel->GetMaximumCS() ) {
	 fChannel->SetMaximumCS(sigma);
       }
       
       if( (sigma < fChannel->GetMinimumCS()) && (sigma > 0.) ) {
	 fChannel->SetMinimumCS(sigma);
       }
       
       ff[ichan] = sigma;
     }
}

//------------------------------------------------------------------------------
//HepMC::GenEvent* Udod::Event()
void Udod::Event()
{
   fChannel = fChannelRegistry->SelectChannel();

   //cout << "Simulate " << fChannel->GetName() << endl;

   if( !fChannel ) {
        cout << " Udod::Event: WARNING: "
                "None of decay channels is set." << endl;
   } 

   while(1) {
      ComputePhaseSpace();
      double sigma = ComputeCrossSection();
      fEventsTotal++;
      double maxCS =  fChannel->GetMaximumCS();
      if( sigma > maxCS ) {
	cout << "New maximum found " << sigma/maxCS << endl;
      }

      if( sigma > Random()*maxCS*1.05 ) { //event accepted
         fEventsAccepted++;

         //-------- Decay and fill decay final state (DFS)
	 fChannel->Decay();

	 UdodHisto::GetInstance()->Fill(1, s1,1);
	 UdodHisto::GetInstance()->Fill(2, s2,1);
	 UdodHisto::GetInstance()->Fill(3, -t1,1);
	 UdodHisto::GetInstance()->Fill(4, -t2,1);

         if( fVerbosity >= DEBUG && fChannel ) {
	   //  cout << "W= " << setprecision(15) << W << endl;
	   UdodHisto::GetInstance()->Fill(6,W,1);
	   //   cout << " Final state: " << endl;
	   for(unsigned int l = 0; l < 3; l++) {
	      // cout << " k# " << l
	      // 	   << " : E_k= " << setprecision(15) << FS[l].Get_E() 
	      // 	   << " P_k= "   << setprecision(15) << FS[l].Get_P() 
	      // 	   << " M_k= "   << setprecision(15) << FS[l].Get_M() 
	      // 	   << " " << FS[l].Theta() << endl;
	      UdodHisto::GetInstance()->Fill(7+l, cos(FS[l].Theta()),1);
	    }

	   //  fChannel->ShowResults();
	    UdodHisto::GetInstance()->Fill(10, cos(DFS[0].Theta()),1);
	    UdodHisto::GetInstance()->Fill(11, cos(DFS[1].Theta()),1);
         }

         // return FillEventRecord();
         return;// 0;
      } else {
         // EffCrossSectCalcBad++;
      }
   } // end of while(1)
}

