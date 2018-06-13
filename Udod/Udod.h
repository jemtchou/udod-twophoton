#ifndef UDOD_h
#define UDOD_h

#include <vector>

//#include "CLHEP/Random/RandomEngine.h"
//#include "HepMC/GenEvent.h"

#include "UdodCut.h"
#include "UdodFinalState.h"
#include "UdodRandomEngine.h"
#include "UdodChannelRegistry.h"

#include "UdodHisto.h"

namespace UDOD
{
  enum {ELECTRON,POSITRON,XSTATE};
  static std::string ParticleName[]= {"e-","e+","X"};

  enum {FATAL,ERROR,WARNING,INFO,DEBUG};
  static std::string DebugLevel[]= {"FATAL","ERROR","WARNING","INFO","DEBUG"};

  enum {VMD, GVMD, RHOPOLE};
  static std::string GammaModel[]= {"VDM","GVDM","Rho-pole"};

  enum {ELECTRON_PAIR, MUON_PAIR, TAU_PAIR};

  class UdodChannelBase;
  class UdodCut;
  class UdodFinalState;

  class Udod
  {
  public:
    Udod(double Ecm);
    ~Udod();

    // set cuts
    void SetElectronEnergyCut(double min, double max) {
      SetEnergyCut(ELECTRON,min,max);
    }

    void SetPositronEnergyCut(double min, double max) {
      SetEnergyCut(POSITRON,min,max);
    }

    void SetElectronThetaCut(double min, double max) {
      SetThetaCut(ELECTRON,min,max);
    }

    void SetPositronThetaCut(double min, double max) {
      SetThetaCut(POSITRON,min,max);
    }

    void SetHadronMassCut(double min, double max);

    void SetVerbosity(int verb) {
      fVerbosity = verb;
    }

    const UdodFinalState* GetFS () const    {
      return FS;
    }
    const std::vector<UdodFinalState>& GetFinalState () const {
      return DFS;
    }

    // calculate total cross-section
    bool Initialize();

    // choose which particles you'd like to product
    //void SetChannel(std::vector<double>& Mass);
    // for couple of particles; num= ELECTRON_PAIR, MUON_PAIR or TAU_PAIR
    // void SetChannel(int num);

    void RegisterChannel(UdodChannelBase* ch);

    const UdodChannelBase* GetChannel() {
      return fChannel;
    }

    // simulate one event
    //       HepMC::GenEvent* Event();
    void Event();

    // finalize
    void Finalize();

  protected:
    void SetEnergyCut(int flag, double min, double max);
    void SetThetaCut(int flag, double min, double max);

    // calculate phase space
    //long double ComputePhaseSpace(const double* rand);
    long double ComputePhaseSpace();

    bool CalculateT1(const double rand);
    bool CalculateT2(const double rand);
    bool CalculateS1(const double rand);
    bool CalculateS2(const double rand);
    bool CalculateMomenta();

    inline long double lambda(long double x,long double y,long double z)
    {
      return (x-y-z)*(x-y-z) - 4*y*z;
    }

    // calculate phasespace limits from kinematic limits
    bool ComputeLimits(double mX, UdodCut& tCut, UdodCut& sCut,
		       UdodCut& thetaCut, UdodCut& energyCut);

    // calculate scattered electron angle
    void ComputeElectronScattering(int particle);
    // calculate properties of two-photon state
    void ComputeGammaGammaState();

    // calculate cross-section
    long double ComputeCrossSection();

    double sigma_ab(int a, int b);
    double ht_or_s_GVMD(double Q2,int transversity);
    double ht_or_s_VMD(double Q2,int transversity);
    double ht_or_s_pole(double Q2,int transversity);

    void Integrand(const double xx[], double* ff);

    double Random();
    void RandomArray(int size, double* rndarray);
    void Message(std::string mess, unsigned int level=1);	 
	 
  private:
    // reset on each run
    int fGammaModel;
    UdodCut fThetaCut[2];
    UdodCut fEnergyCut[2];
    UdodCut fXMassCut;

    double fSigmaMin;
    double fSigmaMax;

    long int fEventsTotal;
    long int fEventsAccepted;

    long double s, rs, beta, beta2; // beam c.m. energy
    
    // reset on each event
    long double W,W2; // two-photon invariant mass
    long double s1,s2,t1,t2;
    long double ds1, ds2, dt1, dt2, dtau;
    UdodCut tCut[2];
    UdodCut sCut[2];
    long double phsp;

    long double E1, E2, EX, P1, P2, PX, sth1, sth2, sthX, cospht, sinpht;

    // final state
    UdodFinalState FS[3];
    std::vector<UdodFinalState> DFS;

    UdodChannelBase* fChannel;
    friend class UdodChannelBase;

    friend class UdodChannel;

    UdodRandomEngine fEngine;
    unsigned int fVerbosity;

    UdodChannelRegistry* fChannelRegistry;
    friend class UdodChannelRegistry;
  };
}

#endif
