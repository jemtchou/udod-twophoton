#ifndef UdodChannelBase_h
#define UdodChannelBase_h

#include <iostream>
#include <vector>

class CkTable;

namespace UDOD
{
  class UdodChannelBase
  {
  public:
    UdodChannelBase(){
      isActive = true;
      maxCS=0;
      integralCS=0;
    }

    virtual ~UdodChannelBase();

    void Initialize(std::vector<double>& Mass);

    double Random(); 

    bool IsActive() { 
      return isActive;
    }

    void Activate() { 
      isActive=true;
    }

    void Inactivate() {
      isActive=false;
    }

    
    int GetID() {
      double ID1,ID2,ID3;
      //ID1 - Leptonic/Hadronic/BW/VMD_BW;
      //ID2 - e+e-/mu+mu-/GVMD/VMDc/rho_pole/rho_pole0;
      //ID3 - Luminosity/Dominick2pi/Dominick2K/BijCorn2pi/Terazawa/pi0;
      
      double ID = ID1*100+ID2*10+ID3;
      //idchan = ID; 
      
      return idchan;

    }

    virtual double Sigma(double s,double s1,double s2,double t1,double t2,double W2,double cphi,double stheta)=0;
    
    void SetIntegralCS(double cs) {integralCS=cs;}
    double GetIntegralCS() {return integralCS;}
    
    void SetMaximumCS(double cs) {maxCS=cs;}
    double GetMaximumCS() {return maxCS;}
    void SetMinimumCS(double cs) {minCS=cs;}
    double GetMinimumCS() {return minCS;}

    unsigned int npart;
    unsigned int GetNparticles(){return npart;}

    void CouplesProduction();
    void KProcedure();
    void Decay() 
    {
      if (npart == 2){CouplesProduction();} 
      else {KProcedure();}
    }

    void ShowResults();
    std::string GetName(){return name;}

  protected:
    std::string name;
       
    std::vector<double> mass_v;
    std::vector<CkTable*> tbls; // tables of integrals

    bool isActive;
    double integralCS;
    double maxCS;
    double minCS;

    int idchan;
  };
};

#endif
