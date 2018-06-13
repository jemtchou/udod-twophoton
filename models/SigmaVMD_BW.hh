#ifndef _SigmaVMD_BW_h
#define _SigmaVMD_BW_h

#include "UdodChannelBase.h"

namespace UDOD{
  
  class SigmaVMD_BW : public UdodChannelBase
  {
  private:
    double mRhoSq;
    double m_M;
    int m_J;
    int m_JP;
    double m_Gamma;
  public:
    SigmaVMD_BW(double M,int J,int JP,double Gamma);
    ~SigmaVMD_BW(){};

    //double GammaGG(int JP);
    double GammaGG();
 
    // BW ------------------------
    //double BW(double W2, double M, double Gamma); 
    double BW(double W2);

    //a,b,Mp - only BW:
    //double sigmaAB_VMD_BW(double W2,double Q12,double Q22,double M,int J,int JP,double Gamma);
    double sigmaAB_VMD_BW(double W2,double Q12,double Q22);
   
    //-----------------------------------------------

    double Sigma(double s,double s1,double s2,double t1,double t2,double W2,double phi,double theta);
  };
}

#endif
