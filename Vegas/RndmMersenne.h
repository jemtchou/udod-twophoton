#ifndef RndmMersenne_h
#define RndmMersenne_h

#include "VegasRandom.h"


class RndmMersenne : public VegasRandom
{/*
	PART 2: Mersenne Twister pseudo-random-long int generator
	adapted from T. Nishimura's and M. Matsumoto's C code at
	http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
*/
 public:
 RndmMersenne(int ndim): VegasRandom(ndim) {};

/* length of state vector */
#define MERSENNE_N 624

/* period parameter */
#define MERSENNE_M 397

#define DEFAULT_SEED 5489

/* 32 or 53 random bits */
#define RANDOM_BITS 32

  void IniRandom(const long int n);
  
  void GetRandom(double *x);
  
  void SkipRandom(long int n);
  
  
 protected:
  unsigned int  MersenneInt(int next);
  
  unsigned int  Twist(unsigned int a, unsigned int b);
  
  void MersenneReload();
  
 private:
  unsigned int m_state[MERSENNE_N];
  int m_next;
};



#endif
