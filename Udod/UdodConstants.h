#ifndef UDOD_CONSTANTS_h
#define UDOD_CONSTANTS_h

  const long double me=0.511E-3; // electron mass, GeV (GALUGA)
  const long double alpha=729927007.0E-11; //QED constant (GALUGA)
//const long double alpha = 7.2973525698e-3;  // QED constant (PDG)
//const long double me   =    0.510998928E-3; // electron mass, GeV (PDG)
//const long double mmu  =  105.6583715E-3;   // muon mass, GeV (PDG)
const long double mtau = 1776.82E-3;        // tau-lepton mass. GeV (PDG)
const long double mmu  =  0.1057; // GALUGA
const double m_rho = 0.7754; // GeV;

const long double me2=me*me;
const long double mmu2=mmu*mmu;
const long double mtau2=mtau*mtau;

const double m_pi=134.9766E-3; // pi^0 mass, GeV (PDG)
// const double m_kaon=493.677E-3; // K^{+/-} mass, GeV (PDG)
// const double f_pi=92.0E-3; // ?? pion decay constant in Gev

const double pi = 3.14159265;// (GALUGA)
			       //const long double pi = 3.1415926535897932384626433832795029L;
const long double to_deg = 180./pi;
const long double to_rad = pi/180;
//const long double nbarn = 0.389379338E6; // GeV-2->nanobarn, PDG
const long double nbarn = 0.389385E+06; // GALUGA
#endif
