#include "RndmMersenne.h"

extern int ndim_;

unsigned int RndmMersenne::Twist(unsigned int a, unsigned int b)
{
  unsigned int mixbits = (a & 0x80000000) | (b & 0x7fffffff);
  unsigned int matrixA = (-(b & 1)) & 0x9908b0df;
  return (mixbits >> 1) ^ matrixA;
}


void RndmMersenne::MersenneReload()
{
  unsigned int *s = m_state;
  int j;

  for( j = MERSENNE_N - MERSENNE_M + 1; --j; ++s )
    *s = s[MERSENNE_M] ^ Twist(s[0], s[1]);
  for( j = MERSENNE_M; --j; ++s )
    *s = s[MERSENNE_M - MERSENNE_N] ^ Twist(s[0], s[1]);
  *s = s[MERSENNE_M - MERSENNE_N] ^ Twist(s[0], m_state[0]);
}


void RndmMersenne::IniRandom(long int lseed)
{
  unsigned int seed = (unsigned int)lseed;
  unsigned int *next = m_state;
  int j;

  for( j = 1; j <= MERSENNE_N; ++j ) {
    *next++ = seed;
    seed = 0x6c078965*(seed ^ (seed >> 30)) + j;
    /* see Knuth TAOCP Vol 2, 3rd Ed, p. 106 for multiplier */
  }

  MersenneReload();
  m_next = 0;
}


unsigned int RndmMersenne::MersenneInt(int next)
{
  unsigned int s = m_state[next];
  s ^= s >> 11;
  s ^= (s << 7) & 0x9d2c5680;
  s ^= (s << 15) & 0xefc60000;
  return s ^ (s >> 18);
}


void RndmMersenne::GetRandom(double *x)
{
  int next = m_next, dim;

  for( dim = 0; dim < ndim_; ++dim ) {
#if RANDOM_BITS == 53
    unsigned int a, b;
#endif

    if( next >= MERSENNE_N ) {
      MersenneReload();
      next = 0;
    }

#if RANDOM_BITS == 53
    a = MersenneInt(next++) >> 5;
    b = MersenneInt(next++) >> 6;
    x[dim] = (67108864.*a + b)/9007199254740992.;
#else
    x[dim] = MersenneInt(next++)/4294967295.;
#endif
  }

  m_next = next;
}

void RndmMersenne::SkipRandom(long int n)
{
#if RANDOM_BITS == 53
  n = 2*n*ndim_ + mersenne_.next;
#else
  n = n*ndim_ + m_next;
#endif
  m_next = n % MERSENNE_N;
  n /= MERSENNE_N;
  while( n-- ) MersenneReload();
}
