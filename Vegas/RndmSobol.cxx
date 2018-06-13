#include "RndmSobol.h"

#include <string.h>

void RndmSobol::IniRandom(const long int n)
{
  static long int ini[9*40] = {
      3,   1,   0,   0,   0,   0,   0,   0,   0,
      7,   1,   1,   0,   0,   0,   0,   0,   0,
     11,   1,   3,   7,   0,   0,   0,   0,   0,
     13,   1,   1,   5,   0,   0,   0,   0,   0,
     19,   1,   3,   1,   1,   0,   0,   0,   0,
     25,   1,   1,   3,   7,   0,   0,   0,   0,
     37,   1,   3,   3,   9,   9,   0,   0,   0,
     59,   1,   3,   7,  13,   3,   0,   0,   0,
     47,   1,   1,   5,  11,  27,   0,   0,   0,
     61,   1,   3,   5,   1,  15,   0,   0,   0,
     55,   1,   1,   7,   3,  29,   0,   0,   0,
     41,   1,   3,   7,   7,  21,   0,   0,   0,
     67,   1,   1,   1,   9,  23,  37,   0,   0,
     97,   1,   3,   3,   5,  19,  33,   0,   0,
     91,   1,   1,   3,  13,  11,   7,   0,   0,
    109,   1,   1,   7,  13,  25,   5,   0,   0,
    103,   1,   3,   5,  11,   7,  11,   0,   0,
    115,   1,   1,   1,   3,  13,  39,   0,   0,
    131,   1,   3,   1,  15,  17,  63,  13,   0,
    193,   1,   1,   5,   5,   1,  27,  33,   0,
    137,   1,   3,   3,   3,  25,  17, 115,   0,
    145,   1,   1,   3,  15,  29,  15,  41,   0,
    143,   1,   3,   1,   7,   3,  23,  79,   0,
    241,   1,   3,   7,   9,  31,  29,  17,   0,
    157,   1,   1,   5,  13,  11,   3,  29,   0,
    185,   1,   3,   1,   9,   5,  21, 119,   0,
    167,   1,   1,   3,   1,  23,  13,  75,   0,
    229,   1,   3,   3,  11,  27,  31,  73,   0,
    171,   1,   1,   7,   7,  19,  25, 105,   0,
    213,   1,   3,   5,   5,  21,   9,   7,   0,
    191,   1,   1,   1,  15,   5,  49,  59,   0,
    253,   1,   1,   1,   1,   1,  33,  65,   0,
    203,   1,   3,   5,  15,  17,  19,  21,   0,
    211,   1,   1,   7,  11,  13,  29,   3,   0,
    239,   1,   3,   7,   5,   7,  11, 113,   0,
    247,   1,   1,   5,   3,  15,  19,  61,   0,
    285,   1,   3,   1,   1,   9,  27,  89,   7,
    369,   1,   1,   3,   7,  31,  15,  45,  23,
    299,   1,   3,   3,   9,   9,  25, 107,  39 };

  int dim, bit, nbits;
  long int max, *pini = ini;

  for( nbits = 0, max = 1; max <= n; max <<= 1 ) ++nbits;
  m_norm = 1./max;

  for( bit = 0; bit < nbits; ++bit )
    m_v[0][bit] = (max >>= 1);

  for( dim = 1; dim < ndim_; ++dim ) {
    long int *pv = m_v[dim], *pvv = pv;
    long int powers = *pini++, j;
    int inibits = -1, bit;
    for( j = powers; j; j >>= 1 ) ++inibits;

    memcpy(pv, pini, inibits*sizeof(*pini));
    pini += 8;

    for( bit = inibits; bit < nbits; ++bit ) {
      long int newv = *pvv, j = powers;
      int b;
      for( b = 0; b < inibits; ++b ) {
        if( j & 1 ) newv ^= pvv[b] << (inibits - b);
        j >>= 1;
      }
      pvv[inibits] = newv;
      ++pvv;
    }

    for( bit = 0; bit < nbits - 1; ++bit )
      pv[bit] <<= nbits - bit - 1;
  }

  m_seq = 0;
  memset(m_prev, 0, (ndim_)*sizeof(*(m_prev)));
}


void RndmSobol::GetRandom(double *x)
{
  long int seq = m_seq++;
  int zerobit = 0, dim;

  while( seq & 1 ) {
    ++zerobit;
    seq >>= 1;
  }
  
  for( dim = 0; dim < ndim_; ++dim ) {
    m_prev[dim] ^= m_v[dim][zerobit];
    x[dim] = m_prev[dim]*m_norm;
  }
}


void RndmSobol::SkipRandom(long int n)
{
  while( n-- ) {
    long int seq = m_seq++;
    int zerobit = 0, dim;

    while( seq & 1 ) {
      ++zerobit;
      seq >>= 1;
    }

    for( dim = 0; dim < ndim_; ++dim )
      m_prev[dim] ^= m_v[dim][zerobit];
  }
}
