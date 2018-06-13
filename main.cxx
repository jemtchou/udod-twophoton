// Header for this module:-
#include <iostream>
#include <vector>
#include <cstdlib>

//#include "CLHEP/Random/Ranlux64Engine.h"
//#include "CLHEP/Random/RandFlat.h"

#include "TTree.h"
#include "TH1F.h"
#include "TFile.h"

#include "Udod.h"
#include "UdodFinalState.h"
#include "UdodConstants.h"

#include "UdodHisto.h"

#include "models/Leptonic.hh"
#include "models/SigmaBW.hh"
#include "models/SigmaVMD_BW.hh"

using namespace std;
using namespace UDOD;

int main()
{
   // Save the random number seeds in the event
  //  CLHEP::HepRandomEngine* engine = new CLHEP::Ranlux64Engine();
  //  const long* seeds = engine->getSeeds();
   vector<int> m_seeds(2);
   // m_seeds[0] = seeds[0];
   // m_seeds[1] = seeds[1];

   UdodHisto* uh = UdodHisto::GetInstance();
   uh->Hist1d(1,"s1",1000,0,100);
   uh->Hist1d(2,"s2",1000,0,100);
   uh->Hist1d(3,"-t1",1000,0,100);
   uh->Hist1d(4,"-t2",1000,0,100);
   uh->Hist1d(5,"sigma",1000,0,1);
   uh->Hist1d(6,"W",1000,0,100);
   uh->Hist1d(7,"e+ctheta",500,-1.2,1.2);
   uh->Hist1d(8,"e-ctheta",500,-1.2,1.2);
   uh->Hist1d(9,"Xctheta",100,-1,1);
   uh->Hist1d(10,"mu+ctheta",100,-1,1);
   uh->Hist1d(11,"mu-ctheta",100,-1,1);
   uh->Hist1d(20,"s1min",1000,0,100);
   uh->Hist1d(21,"s1max",1000,0,100);
   uh->Hist1d(22,"s2min",1000,0,100);
   uh->Hist1d(23,"s2max",1000,0,100);


   double Ecm = 200; // for LEP

   Udod gena(Ecm);

   gena.SetVerbosity(DEBUG);//output verbosity by default
       
   gena.SetElectronEnergyCut(0,200);
   gena.SetPositronEnergyCut(0,200);
   gena.SetElectronThetaCut(0,180);
   gena.SetPositronThetaCut(0,180);
   gena.SetHadronMassCut(0.22,20);

   gena.RegisterChannel(new Leptonic(mmu2));
   gena.RegisterChannel(new Leptonic(me2));
   gena.RegisterChannel(new SigmaBW(0.548,0.775,0,2,0.00000131));
   gena.RegisterChannel(new SigmaVMD_BW(0.548,0,2,0.00000131));
 
   gena.Initialize();

   int iev = 0;
//   for(iev = 0; iev < 100000; iev++) {
//     cout << " Event# " << iev << endl;
     gena.Event();
//   }

   // ----------- test for pair decay
   // test_pairs(gena);

   // ----------- test for K-procedure
   //Set masses in GeV
   // vector<double> Masses(3,me);
   // gena.SetChannel(Masses);

   gena.Finalize();

   uh->Save();

   //  delete engine;
   return 0;
};

