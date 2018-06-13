#ifndef UdodRandomEngine_h
#define UdodRandomEngine_h

#include "TRandom3.h"

#include <iostream>

namespace UDOD
{
  class UdodRandomEngine
  {
  public:
    UdodRandomEngine() {
      engine = new TRandom3;
    }

    virtual ~UdodRandomEngine() {
      delete engine;
    }

    void SetSeed(int seed) {
      engine->SetSeed(seed);
    }
	 
    int GetSeed() {
      return engine->GetSeed();
    }

    double Random() {
      return engine->Rndm();
    }

    double Random(double a, double b) {
      return a+(b-a)*Random();
    }

    void RandomArray(int n, double* arr) {
      engine->RndmArray(n,arr);
    }

    // Access methods
  private:
    TRandom* engine;
  };
}

#endif
