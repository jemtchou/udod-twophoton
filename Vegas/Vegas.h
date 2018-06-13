#ifndef Vegas_h
#define Vegas_h

//
//  Vegas Monte-Carlo integration
//
//  Based on the Cuba library by Thomas Hahn 
//  Comput. Phys. Commun. 168 (2005) 78 [hep-ph/0404043])
//  Cuba is available at http://www.feynarts.de/cuba/)
//
//  C++ version prepared by A.Zhemchugov (zhemchugov@jinr.ru)
//  Last modified 27 Nov 2008
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cmath>
#include <unistd.h>
#include <fcntl.h>

#include "VegasUtils.h"

class VegasRandom;

// demo defines

#define NDIM 5 
#define NCOMP 1
#define NBINS 128

#define LAST 4

#define MAXGRIDS 10
#define MAXSTATESIZE 128

typedef double Grid[NBINS];

struct Cumulants{
  double sum, sqsum;
  double weightsum, avgsum;
  double chisum, chisqsum, guess;
  double avg, err, chisq;
};

//void Integrand(const int *ndim, const double xx[],
//	       const int *ncomp, double ff[],const double*);

class Integrable
{
 public:
  virtual void Integrand(const int *ndim, const double xx[],
			const int *ncomp, double ff[])=0;
};  


struct State{
  int niter;
  long int nsamples, neval;
  Cumulants cumul[NCOMP];
  Grid grid[NDIM];
};

// Integrate over the unit hypercube

enum {MERSENNE, SOBOL};

class Vegas
{
 public:
  Vegas();
  Vegas(int ndim, int ncomp, int pseudorng);
  ~Vegas();

  void SetFunction(Integrable* func){m_function = func;};

  int Integrate(double *integral, double *error, double *prob, long int *pneval);

 protected:
  void DoSample(long int n, const double *w, const double *x, double *f);
  
  bool BadDimension(const int ndim);
  
  bool BadComponent(const int ncomp);
  
  void GetGrid(Grid *grid);

  void PutGrid(Grid *grid);
  
  void RefineGrid(Grid grid, Grid margsum);

  double MaxErr(double avg);
  
 private:
  VegasRandom* m_rndm;

  Integrable* m_function;

  int m_ndim;
  int m_ncomp;
  long int m_neval;

  double m_epsrel, m_epsabs;

  long int m_mineval, m_maxeval, m_nstart, m_nincrease;
  int m_verbose;

  Grid *gridptr_[MAXGRIDS];
  int griddim_[MAXGRIDS];
  int vegasnbatch;
  int vegasgridno;
  char vegasstate[MAXSTATESIZE];

  int SHARPEDGES;
  int REGIONS;

};



#endif
