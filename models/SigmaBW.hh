#ifndef _SigmaBW_h
#define _SigmaBW_h

#include "UdodChannelBase.h"

namespace UDOD{
  
  class SigmaBW : public UdodChannelBase
  {
  private:
    double m_M;
    double m_Mp;
    int m_J;
    int m_JP;
    double m_Gamma;
    
  public:
    SigmaBW(double M,double Mp,int J,int JP,double Gamma) {
      idchan = 2;
      m_M = M;
      m_Mp = Mp;
      m_J = J;
      m_JP = JP;
      m_Gamma = Gamma;
    };
    ~SigmaBW(){};
  
    // GammaGG
    double GammaGG();

    double nu(double W2,double t1,double t2);
  
    // fAB-----------------------------
    double X(double W2,double Q12,double Q22);
    double kappa(double W2,double Q12,double Q22);
    //double kappa(double W2,double Q12,double Q22, double M2);
  
    // f_TT, J=0-, 0+, 1+, 2+, 2-
    //double fTT_0M(double W2,double Q12,double Q22, double M2);
    double fTT_0M(double W2,double Q12,double Q22);
    double fTT_0P(double W2,double Q12,double Q22);
    double fTT_1P(double W2,double Q12,double Q22);
    double fTT_2P(double W2,double Q12,double Q22);
    double fTT_2M(double W2,double Q12,double Q22);

    // f_TS, J=0-, 0+, 1+, 2+, 2-
    double fTS_0M(double W2,double Q12,double Q22);
    double fTS_0P(double W2,double Q12,double Q22);
    double fTS_1P(double W2,double Q12,double Q22);
    double fTS_2P(double W2,double Q12,double Q22);
    double fTS_2M(double W2,double Q12,double Q22);
  
    // f_ST, J=0-, 0+, 1+, 2+, 2-
    double fST_0M(double W2,double Q12,double Q22);
    double fST_0P(double W2,double Q12,double Q22);
    double fST_1P(double W2,double Q12,double Q22);
    double fST_2P(double W2,double Q12,double Q22);
    double fST_2M(double W2,double Q12,double Q22);

    // f_SS, J=0-, 0+, 1+, 2+, 2-
  
    double fSS_0M(double W2,double Q12,double Q22);
    double fSS_0P(double W2,double Q12,double Q22);
    double fSS_1P(double W2,double Q12,double Q22);
    double fSS_2P(double W2,double Q12,double Q22);
    double fSS_2M(double W2,double Q12,double Q22);
    //------------------------------------

    //-------------------fAB------------------
    //double fAB(int a, int b, double W2,double Q12,double Q22, double Mp, int JP);
    double fAB(int a, int b, double W2,double Q12,double Q22);

    // BW ------------------------
    //double BW(double W2, double M, double Gamma);
    double BW(double W2);
    

    //a,b,Mp - only BW:
    //double sigmaAB_BW(int a, int b,double W2,double Q12,double Q22, double M, double Mp, int J, int JP,double Gamma);
     double sigmaAB_BW(int a, int b,double W2,double Q12,double Q22);
    //-----------------------------------------------

    double Sigma(double s,double s1,double s2,double t1,double t2,double W2,double phi,double theta);
    
  };

}

#endif
