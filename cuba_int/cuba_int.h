#ifndef CUBA_INT_H
#define CUBA_INT_H

#include <string>
#include <vector>

#include "cuba.h"
#include "FastDelegate.h"

enum INTMETH {VEGAS = 1, SUAVE, DIVONNE, CUHRE};

// function with one parameter and return one double value
typedef fastdelegate::FastDelegate2<const double*, double*, void> Integ_Fun;

class cuba_int
{
   public:
      cuba_int();

      // Call integration
      void Integral();
      void Integral(INTMETH nom1)   {
         set_nom(nom1);
         Integral();
      }

      // Set/get integrands
      void set_fun(Integ_Fun func) {
         fun_delegate.clear();
         fun_delegate.push_back(func);
      }
//      void add_fun(Integ_Fun func)  {
//         fun_delegate.push_back(func);
//      }
      const Integ_Fun& get_fun() const {
         return fun_delegate[0];
      }
      // the number of components of the integrand
      unsigned int get_ncomp() const {
	   return ncomp;
//         return fun_delegate.size();
      }
      void set_ncomp(int nc) {ncomp = nc;}

      // access functions
      int get_verbose() const {
         return verbose;
      }
      void set_verbose(int ve) {
         verbose = ve;
      }
      unsigned int get_flags() const {
         return flags;
      }
      void set_flags(unsigned int bit0, unsigned int bit2,
                     unsigned int bit4, unsigned int bit8) {
         flags |= bit2 << 2;
         flags |= bit4 << 4;
         flags |= bit8 << 8;
      }

      // Common parameters
      INTMETH get_nom() const {
         return nom;
      }
      void set_nom(INTMETH nom1) {
         nom = nom1;
      }
      const std::string& str_nom() const {
         static std::string name[] = {"","VEGAS","SUAVE","DIVONNE","CUHRE"};
         return name[nom];
      }
      int get_dim() const {
         return ndim;
      }
      void set_dim(int dim) {
         ndim = dim;
      }
      double get_epsrel() const {
         return epsrel;
      }
      void set_epsrel(double el) {
         epsrel = el;
      }
      double get_epsabs() const {
         return epsabs;
      }
      void set_epsabs(double es) {
         epsabs = es;
      }
      int get_seed() const {
         return seed;
      }
      void set_seed(int sd) {
         seed = sd;
      }
      int get_mineval() const {
         return mineval;
      }
      void set_mineval(int mil) {
         mineval = mil;
      }
      int get_maxeval() const {
         return maxeval;
      }
      void set_maxeval(int mal) {
         maxeval = mal;
      }

      // Vegas-specific parameters
      int get_nstart() const {
         return nstart;
      }
      void set_nstart(int nt) {
         nstart = nt;
      }
      int get_nincrease() const {
         return nincrease;
      }
      void set_nincrease(int ne) {
         nincrease = ne;
      }
      int get_nbatch() const {
         return nbatch;
      }
      void set_nbatch(int nc) {
         nbatch = nc;
      }
      int get_gridno() const {
         return gridno;
      }
      void set_gridno(int gno) {
         gridno = gno;
      }
      const std::string& get_statefile() const {
         return statefile;
      }
      void set_statefile(const std::string& se) {
         statefile = se;
      }

      // Suave-specific parameters
      int get_nnew() const {
         return nnew;
      }
      void set_nnew(int nw) {
         nnew = nw;
      }
      double get_flatness() const {
         return flatness;
      }
      void set_flatness(double fs) {
         flatness = fs;
      }

      // Divonne-specific parameters
      int get_key1() const {
         return key1;
      }
      void set_key1(int k1) {
         key1 = k1;
      }
      int get_key2() const {
         return key2;
      }
      void set_key2(int k2) {
         key2 = k2;
      }
      int get_key3() const {
         return key3;
      }
      void set_key3(int k3) {
         key3 = k3;
      }
      int get_maxpass() const {
         return maxpass;
      }
      void set_maxpass(int ms) {
         maxpass = ms;
      }
      double get_border() const {
         return border;
      }
      void set_border(double br) {
         border = br;
      }
      double get_maxchisq() const {
         return maxchisq;
      }
      void set_maxchisq(double mq) {
         maxchisq = mq;
      }
      double get_mindeviation() const {
         return mindeviation;
      }
      void set_mindeviation(double mn) {
         mindeviation = mn;
      }
      int get_ngiven() const {
         return ngiven;
      }
      void set_ngiven(int nn) {
         ngiven = nn;
      }
      int get_ldxgiven() const {
         return ldxgiven;
      }
      void set_ldxgiven(int ln) {
         ldxgiven = ln;
      }
      void set_xgiven(double xn[]) {
         xgiven = xn;
      }
      void set_peakfinder(peakfinder_t pf) {
         peakfinder = pf;
      }
      int get_nextra() const {
         return nextra;
      }
      void set_nextra(int na) {
         nextra = na;
      }

      // Cuhre-specific parameters
      int get_key() const {
         return key;
      }
      void set_key(int k) {
         key = k;
      }

      // Output variables
      int get_nregions() const {
         return nregions;
      }
      int get_neval() const {
         return neval;
      }
      int get_fail() const {
         return fail;
      }
      double get_integral(unsigned int i=0) const {
         return integral[i];
      }
      double get_error(unsigned int i=0) const {
         return error[i];
      }
      double get_prob(unsigned int i=0) const {
         return prob[i];
      }

   private:
      int verbose; // verbosity level
      int flags;   // flags governing the integration (see manual)

      // Integrands
      std::vector<Integ_Fun> fun_delegate;

      // Common arguments
      INTMETH nom; // integration method
      int ndim;    // dimension of space
      unsigned int ncomp;   // number of components

      double epsrel; // the requested relative accuracy
      double epsabs; // the requested absolute accuracy

      int mineval; // the minimum number of integrand evaluations required
      int maxeval; // the maximum number of integrand evaluations allowed

      int seed; // the seed for the pseudo-random-number generator

      // Vegas-specific arguments
      int nstart; // the initial number of integrand evaluations per iteration
      int nincrease; /* the increase in the number of
                        integrand evaluations per iteration */
      int nbatch; // the batch size for sampling
      int gridno; // the slot in the internal grid table (do not use it)

      std::string statefile; // a filename for storing the internal state

      // Suave-specific arguments
      int nnew; // the number of new integrand evaluations in each subdivision
      double flatness; /* the parameter used to compute the fluctuation
                          of a sample (do not use too large value) */

      // Divonne-specific arguments
      int key1;   // determines sampling in the partitioning phase
      int key2;   // determines sampling in the final integration phase
      int key3;   // sets the strategy for the refinement phase
      int maxpass;   // controls the thoroughness of the partitioning phase
      double border; // the width of the border of the integration region

      double maxchisq; /* the maximum chi^2 value a single subregion is
                          allowed to have in the final integration phase */
      double mindeviation; /* a bound, given as the fraction of the
                              requested error of the entire integral,
                              which determines whether it is worthwhile
                              further examining a region that failed
                              the chi^2 test. */

      int ngiven;     // the number of points in the xgiven array
      int ldxgiven;   // the leading dimension of xgiven
      double* xgiven; // list of points where the integrand might have peaks

      int nextra; /* the maximum number of extra points the peak-finder
                     function will return */
      peakfinder_t peakfinder;   // the peak-finder function

      // Cuhre-specific arguments
      int  key;   // chooses the basic integration rule

      //Output variables
      int nregions;  // the actual number of subregions needed
      int neval;     // the actual number of integrand evaluations
      int fail;      // an error flag: (0 - the desired accuracy was reached)

      std::vector<double> integral; // the integral
      std::vector<double> error;    // the presumed absolute error of integral
      std::vector<double> prob;     /* the chi^2 probability that error
                                       is not a reliable estimate
                                       of the true integration error */
};
#endif
