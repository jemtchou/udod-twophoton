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


#include "Vegas.h"
#include "RndmSobol.h"
#include "RndmMersenne.h"

#include <iostream>
#include <sys/stat.h> // for stat()

Vegas::Vegas(int Ndim, int Ncomp, int Pseudorng)
{
  m_ndim = Ndim;
  m_ncomp = Ncomp;

  if( Pseudorng == SOBOL )
    m_rndm = new RndmSobol(m_ndim);
  else
    m_rndm = new RndmMersenne(m_ndim);

  m_epsrel = 1e-3;
  m_epsabs = 1e-12;
  m_mineval = 0;
  m_maxeval = 1000000;
  m_nstart = 10000;
  m_nincrease = 1000;

  m_verbose = 0;
  
  vegasnbatch = 1000;
  vegasgridno = 0;
  sprintf(vegasstate,"");

  SHARPEDGES=0;
  REGIONS=0;
}

Vegas::~Vegas()
{
  delete m_rndm;
}

double Vegas::MaxErr(double avg)
{ 
  return VegasUtils::Max(m_epsrel*fabs(avg), m_epsabs);
};

void  Vegas::DoSample(long int n, const double *w, const double *x, double *f)
{
  m_neval += n;
  while( n-- ) {
    m_function->Integrand(&m_ndim, x, &m_ncomp, f);//, w++);
    x += m_ndim;
    f += m_ncomp;
  }
}

/* Note: unsigned char must be wide enough to hold the long ints 0..NBINS */

bool Vegas::BadDimension(const int ndim)
{
  //#if NDIM > 0
  //if( m_ndim > NDIM ) return true;
  //#endif
  //return m_ndim < SOBOL_MINDIM || (!m_pseudorng && m_ndim > SOBOL_MAXDIM);
  return false;
}


bool Vegas::BadComponent(const int ncomp)
{
  if (m_ncomp > 0){
    if( ncomp > m_ncomp ) return true;
  }
  else
    return ncomp < 1;

  return false;
}

void Vegas::GetGrid(Grid *grid)
{
  int bin, dim;
  unsigned const int slot = vegasgridno - 1;

  if( slot < MAXGRIDS && gridptr_[slot] ) {
    if( griddim_[slot] == m_ndim ) {
      memcpy(grid, gridptr_[slot], (m_ndim)*sizeof(*grid));
      return;
    }
    free(gridptr_[slot]);
    gridptr_[slot] = NULL;
  }

  for( bin = 0; bin < NBINS; ++bin )
    grid[0][bin] = (bin + 1)/(double)NBINS;
  for( dim = 1; dim < m_ndim; ++dim )
    memcpy(&grid[dim], &grid[0], sizeof(grid[dim]));
}

/*********************************************************************/

void Vegas::PutGrid(Grid *grid)
{
  unsigned const int slot = vegasgridno - 1;

  if( slot < MAXGRIDS ) {
    if( gridptr_[slot] == NULL ) 
      gridptr_[slot] = (Grid*)malloc(m_ndim*sizeof(gridptr_[slot]));
    if( gridptr_[slot] == NULL )
      {
	std::cerr << "PutGrid: Out of memory" << std::endl;
	exit(1);
      }
    griddim_[slot] = m_ndim;
    memcpy(gridptr_[slot], grid, (m_ndim)*sizeof(*(gridptr_[slot])));
  }
}

/*********************************************************************/

void Vegas::RefineGrid(Grid grid, Grid margsum)
{
  double avgperbin, thisbin, newcur, delta;
  Grid imp, newgrid;
  int bin, newbin;

  /* smooth the f^2 value stored for each bin */
  double prev = margsum[0];
  double cur = margsum[1];
  double norm = margsum[0] = .5*(prev + cur);
  for( bin = 1; bin < NBINS - 1; ++bin ) {
    const double s = prev + cur;
    prev = cur;
    cur = margsum[bin + 1];
    norm += margsum[bin] = (s + cur)/3.;
  }
  norm += margsum[NBINS - 1] = .5*(prev + cur);

  if( norm == 0 ) return;
  norm = 1/norm;

  /* compute the importance function for each bin */
  avgperbin = 0;
  for( bin = 0; bin < NBINS; ++bin ) {
    double impfun = 0;
    if( margsum[bin] > 0 ) {
      const double r = margsum[bin]*norm;
      avgperbin += impfun = pow((r - 1)/log(r), 1.5);
    }
    imp[bin] = impfun;
  }
  avgperbin /= NBINS;

  /* redefine the size of each bin */
  cur = newcur = 0;
  thisbin = 0;
  bin = -1;
  for( newbin = 0; newbin < NBINS - 1; ++newbin ) {
    while( thisbin < avgperbin ) {
      thisbin += imp[++bin];
      prev = cur;
      cur = grid[bin];
    }
    thisbin -= avgperbin;
    delta = (cur - prev)*thisbin;
    newgrid[newbin] = SHARPEDGES ?
      cur - delta/imp[bin] :
      (newcur = VegasUtils::Max(newcur,
        cur - 2*delta/(imp[bin] + imp[VegasUtils::IDim(bin - 1)])));
  }
  memcpy(grid, newgrid, (NBINS-1)*sizeof(*(grid)));

  grid[NBINS - 1] = 1;
}

int Vegas::Integrate(double *integral, double *error, double *prob, long int *pneval)
{
  double *sample;
  int dim, comp;
  int fail = 1;

  int statemsg = m_verbose;
  struct stat st;

  State state;

  if( BadComponent(m_ncomp) || BadDimension(m_ndim) ) 
    return -1;
    
  m_neval = 0;
  
  if( m_verbose > 0 ) {
    std::cout << "Vegas input parameters:\n"
	      << "  Number of dimensions of the integral\t" << m_ndim << std::endl 
	      << "  Number of components of the integrand\t" << m_ncomp << std::endl
	      << "  Relative accuracy " << m_epsrel << std::endl 
	      << "  Absolute accuracy " << m_epsabs << std::endl
	      << "  Minimal number of integrand evaluations\t" << m_mineval << std::endl 
	      << "  Maximal number of integrand evaluations\t" << m_maxeval << std::endl
	      << "  Initial number of evaluations per iteration\t" << m_nstart << std::endl 
	      << "  Increase of the number of evaluations per iteration\t" << m_nincrease << std::endl
	      << "  vegasgridno " << vegasgridno << std::endl 
	      << "  vegasstate \"" << vegasstate << "\"" << std::endl;
  }
  
  m_rndm->IniRandom(2*m_maxeval);

  if( *vegasstate && stat(vegasstate, &st) == 0 &&
      st.st_size == sizeof(state) && (st.st_mode & 0400) ) {
    const int h = open(vegasstate, O_RDONLY);
    read(h, &state, sizeof(state));
    close(h);
    m_rndm->SkipRandom(m_neval = state.neval);

    if( m_verbose ) 
      std::cout << "Restoring state from " << vegasstate << std::endl;
  }
  else {
    state.niter = 0;
    state.nsamples = m_nstart;
    memset(state.cumul, 0, sizeof(state.cumul));
    GetGrid(state.grid);
  }

  sample = (double*)malloc((vegasnbatch)*((m_ndim + m_ncomp + 1)*sizeof(double) + m_ndim*sizeof(unsigned char)));
 
  if(sample == 0)
    std::cerr << "Integrate: Malloc failed " << std::endl;
    
  /* main iteration loop */

  for( ; ; ) {
    long int nsamples = state.nsamples;
    const double jacobian = 1./nsamples;
    Grid margsum[NCOMP][NDIM];
    
    memset(margsum, 0, sizeof(margsum));
   
    for( ; nsamples > 0; nsamples -= vegasnbatch ) {
      const long int nbatch = VegasUtils::IMin(vegasnbatch, nsamples);
      double *w = sample;
      double *x = w + nbatch;
      double *f = x + nbatch*m_ndim;
      double *lastf = f + nbatch*m_ncomp;
      unsigned char *bin = (unsigned char *)lastf;
      
      while( x < f ) {
        double weight = jacobian;

        m_rndm->GetRandom(x);
	for( dim = 0; dim < m_ndim; ++dim ) {
          const double pos = *x*NBINS;
          const int ipos = (int)pos;
          const double prev = (ipos == 0) ? 0 : state.grid[dim][ipos - 1];
          const double diff = state.grid[dim][ipos] - prev; 
          *x++ = prev + (pos - ipos)*diff;
	  *bin++ = ipos;
	  weight *= diff*NBINS;
        }

        *w++ = weight;
      }

      DoSample(nbatch, sample, w, f);

      w = sample;
      bin = (unsigned char *)lastf;

      while( f < lastf ) {
        const double weight = *w++;

        for( comp = 0; comp < m_ncomp; ++comp ) {
          double wfun = weight*(*f++);
          if( wfun ) {
            Cumulants *c = &state.cumul[comp];
            Grid *m = margsum[comp];

            c->sum += wfun;
            c->sqsum += wfun *= wfun;
            for( dim = 0; dim < m_ndim; ++dim )
              m[dim][bin[dim]] += wfun;
          }
        }

        bin += m_ndim;
      }
    }

    fail = 0;

    /* compute the integral and error values */
    
    for( comp = 0; comp < m_ncomp; ++comp ) {
      Cumulants *c = &state.cumul[comp];
      double avg, sigsq;
      double w = VegasUtils::Weight(c->sum, c->sqsum, state.nsamples);
      
      sigsq = 1/(c->weightsum += w);
      avg = sigsq*(c->avgsum += w*c->sum);
      
      c->avg = LAST ? (sigsq = 1/w, c->sum) : avg;
      c->err = sqrt(sigsq);
      fail |= (c->err > MaxErr(c->avg));

      if( state.niter == 0 ) c->guess = c->sum;
      else {
        c->chisum += w *= c->sum - c->guess;
        c->chisqsum += w*c->sum;
      }
      c->chisq = c->chisqsum - avg*c->chisum;

      c->sum = c->sqsum = 0;
    }

    if( m_verbose ) {
      std::cout << "Iteration " << (state.niter+1) << ":  " 
		<< m_neval << " integrand evaluations so far" << std::endl;
      
      for( comp = 0; comp < m_ncomp; ++comp ) {
        const Cumulants *c = &state.cumul[comp];
        std::cout << "[" << (comp+1) << "] " << c->avg << " +- " 
	     << c->err << "  \tchisq " << c->chisq << " (" 
	     << state.niter << " df)" << std::endl;
      }
    }

    if( fail == 0 && m_neval >= m_mineval ) {
      if( *vegasstate ) unlink(vegasstate);
      break;
    }

    if( m_neval >= m_maxeval && *vegasstate == 0 ) break;

    if( m_ncomp == 1 )
      for( dim = 0; dim < m_ndim; ++dim )
        RefineGrid(state.grid[dim], margsum[0][dim]);
    else {
      for( dim = 0; dim < m_ndim; ++dim ) {
        Grid wmargsum;
	memset(wmargsum, 0, sizeof(wmargsum));
	for( comp = 0; comp < m_ncomp; ++comp ) {
          double w = state.cumul[comp].avg;
          if( w != 0 ) {
            const double *m = margsum[comp][dim];
            int bin;
            w = 1/ VegasUtils::Sq(w);
            for( bin = 0; bin < NBINS; ++bin )
              wmargsum[bin] += w*m[bin];
          }
        }
        RefineGrid(state.grid[dim], wmargsum);
      }
    }

    ++state.niter;
    state.nsamples += m_nincrease;

    if( *vegasstate ) {
      const int h = creat(vegasstate, 0666);
      if( h != -1 ) {
        state.neval = m_neval;
        write(h, &state, sizeof(state));
        close(h);

        if( statemsg ) {
	  std::cout << "\nSaving state to " << vegasstate << std::endl;
	  statemsg = 0;
        }
      }
      if( m_neval >= m_maxeval ) break;
    }
  }

  for( comp = 0; comp < m_ncomp; ++comp ) {
    const Cumulants *c = &state.cumul[comp];
    integral[comp] = c->avg;
    error[comp] = c->err;
    prob[comp] =  VegasUtils::ChiSquare(c->chisq, state.niter);
  }

  free(sample);
  PutGrid(state.grid);

  *pneval = m_neval;

  return fail;
};
