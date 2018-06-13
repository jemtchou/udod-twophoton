#ifndef RndmSobol_h
#define RndmSobol_h

#include "VegasRandom.h"

#define SOBOL_MINDIM 1
#define SOBOL_MAXDIM 40

class RndmSobol : public VegasRandom
{
  /*
    PART 1: Sobol quasi-random-number generator
    adapted from ACM TOMS algorithm 659
  */
 public:
 RndmSobol(int ndim): VegasRandom(ndim) {};
  void IniRandom(const long int n);
  
  void GetRandom(double *x);
  
  void SkipRandom(long int n);

 private:
  double m_norm;
  long int m_v[SOBOL_MAXDIM][30], m_prev[SOBOL_MAXDIM];
  long int m_seq;  
};

#endif
