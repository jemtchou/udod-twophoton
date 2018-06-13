#include <iostream>
#include <cstdlib>
#include <cassert>

#include "cuba_int.h"

using namespace std;
using namespace fastdelegate;

extern "C" {
   // auxiliary function for interface to Cuba library
   //*************************************************************************80
   static int IntegFD(const int* ndim, const double xx[],
                      const int* ncomp, double ff[], void* userdata)
   //*************************************************************************80
   {
      const cuba_int* cub = reinterpret_cast<const cuba_int*>(userdata);
      assert( cub );

//      for(int icomp = 0; icomp < *ncomp; icomp++) {
//         ff[icomp] = (cub->get_fun(icomp))(xx);
//      }
      (cub->get_fun())(xx,ff);

      return 0;
   }
} // end of extern "C"

//****************************************************************************80
cuba_int::cuba_int()
//****************************************************************************80
{
   // Common arguments
   nom = VEGAS;
   ndim = 5;
   verbose = 0;
   flags = 0;

   // the requested relative and absolute accuracies:
   epsrel = 1e-3;
   epsabs = 1e-12;

   seed = 0; // Use Sobol (quasi-random) generator

   // the minimum and maximum numbers of integrand evaluations:
   mineval = 0;
   maxeval = 1000000;

   // Vegas-specific arguments
   nstart = 1000;
   nincrease = 500;
   nbatch = 1000; /* nbatch = 1000 is a reasonable value, though
                     it should not affect performance too much */
   gridno = 0;    // Do not use it!
   // statefile.clean();

   // Suave-specific arguments
   nnew = 1000;
   flatness = 25.;

   // Divonne-specific arguments
   key1 = 47;
   key2 = 1;
   key3 = 1;
   maxpass = 5;
   border = 0.;
   maxchisq = 10.;
   mindeviation = 0.25;

   ngiven = 0;
   ldxgiven = ndim;
   xgiven = 0;

   nextra = 0;
   peakfinder = 0;

   // Cuhre-specific arguments
   key = 0; /* use the default key=0: the degree-13 rule in 2 dimensions,
               the degree-11 rule in 3 dimensions,
               and the degree-9 rule otherwise */
}

//****************************************************************************80
void cuba_int::Integral()
//****************************************************************************80
{
   int gridno = 0;   // the slot in the internal grid table. Do not use it.
   void* userdata = this;  // user data passed to the integrand.

   // the number of components of the integrand
   //unsigned int ncomp = get_ncomp();
   cout << "NCOMP " << ncomp << endl;
   integral.clear();
   integral.resize(ncomp,0);
   error.clear();
   error.resize(ncomp,0);
   prob.clear();
   prob.resize(ncomp,0);

   neval = 0;
   fail = 0;

   switch( nom ) {
   case VEGAS:
      Vegas(   ndim, ncomp, IntegFD, userdata, epsrel, epsabs,
               flags, seed, mineval, maxeval, nstart, nincrease,
               nbatch, gridno, statefile.c_str(), &neval, &fail,
               &integral[0], &error[0], &prob[0]
           );
      if( verbose ) {
         cout << " vegas:\tneval " << neval << "\tfail " << fail << endl;

         for(unsigned int icomp = 0; icomp < ncomp; icomp++) {
            cout << " vegas:\t" << integral[icomp] << " +/- "
                 << error[icomp] << "\t prob= " << prob[icomp] << endl;
         }
      }
      break;

   case SUAVE:
      Suave(   ndim, ncomp, IntegFD, userdata, epsrel, epsabs,
               flags, seed, mineval, maxeval, nnew, flatness, &nregions,
               &neval, &fail, &integral[0], &error[0], &prob[0]
           );
      if( verbose ) {
         cout << " suave:\tnregions " << nregions
              << "\tneval " << neval << "\tfail " << fail << endl;

         for(unsigned int icomp = 0; icomp < ncomp; icomp++) {
            cout << " suave:\t" << integral[icomp] << " +/- "
                 << error[icomp] << "\t prob= " << prob[icomp] << endl;
         }
      }
      break;

   case DIVONNE:
      Divonne( ndim, ncomp, IntegFD, userdata, epsrel, epsabs, flags,
               seed, mineval, maxeval, key1, key2, key3, maxpass,
               border, maxchisq, mindeviation, ngiven, ldxgiven, xgiven,
               nextra, peakfinder, &nregions, &neval, &fail,
               &integral[0], &error[0], &prob[0]
             );
      if( verbose ) {
         cout << " divonne: nregions " << nregions
              << "\tneval " << neval << "\tfail " << fail << endl;

         for(unsigned int icomp = 0; icomp < ncomp; icomp++) {
            cout << " divonne: " << integral[icomp] << " +/- "
                 << error[icomp] << "\t prob= " << prob[icomp] << endl;
         }
      }
      break;

   case CUHRE:
      Cuhre(  ndim, ncomp, IntegFD, userdata,   epsrel, epsabs, flags,
              mineval, maxeval, key, &nregions, &neval, &fail,
              &integral[0], &error[0], &prob[0]
           );
      if( verbose ) {
         cout << " cuhre:\tnregions " << nregions
              << "\tneval " << neval << "\tfail " << fail << endl;

         for(unsigned int icomp = 0; icomp < ncomp; icomp++) {
            cout << " cuhre:\t" << integral[icomp] << " +/- "
                 << error[icomp] << "\t prob= " << prob[icomp] << endl;
         }
      }
      break;

   default:
      cout << "cuba_int::Integral: ERROR! WRONG 'nom'" << endl;
      exit(1);
   }
}

#ifdef MAINCUBAINT
#include <numeric>
#include <cmath>
#include <cstdio>

//****************************************************************************80
class testfun
{
   private:
      int ndim; // dimension of space
      double beta; // parameter of Oscillatory fun.
      double max_fun; // save here max of fun.

   public:
      testfun(int n, double b) : ndim(n), beta(b) {
         max_fun=0;
      }

      double MaxFun() {
         return max_fun;
      }
      void ResetMaxFun() {
         max_fun = 0.;
      }

      double Oscillatory_F(const double z[]) {
         double sum   = accumulate(z, z+ndim, 0.);
         double total = 2*M_PI*beta + sum;
         double value = cos(total);

         if( value > max_fun ) {
            max_fun = value;
            // cout << " new max " << max_fun << endl;
         }

         return value;
      }

      double Oscillatory_I() {
//
// \integral_{0}^{1} [cos(2*pi*beta + \sum_{i=1}^{Ndim} X_i] d\vec{X}=
// (2*sin(0.5))^Ndim * cos(2*pi*beta + Ndim/2)
//
         static const double tmp = 2*sin(0.5);
         double value = pow(tmp,ndim) * cos(2.0*M_PI*beta + 0.5*ndim);

         return value;
      }

      void InternalIntegrator();
};
//****************************************************************************80

//****************************************************************************80
void testfun::InternalIntegrator()
//****************************************************************************80
{
   // the exact value of the integral:
   double ex_val = Oscillatory_I();
   // cout << " accurate integral= " << ex_val << endl;

   // numeric integration:
   cuba_int cuba;
   cuba.set_dim(ndim);
   cuba.set_fun( MakeDelegate(this, &testfun::Oscillatory_F) );

   // try different integration methods
   printf("   type      value        abs.err"
          "      true val    true rel.err"
          "   max.fun\n");
   for(int j = VEGAS; j <= CUHRE; j++) { // VEGAS = 1, SUAVE, DIVONNE, CUHRE
      cuba.set_nom( (INTMETH)j );
      cuba.Integral();
      double val = cuba.get_integral();
      double abserr = cuba.get_error();
      double eps = fabs(val/ex_val - 1.);

      printf(" %7s %12.6f +/- %8.6f %12.6f %12.2e  %12.6f\n",
             (cuba.str_nom()).c_str(),val,abserr,ex_val,eps,MaxFun());

      ResetMaxFun();
   }

}

//****************************************************************************80
int main()
//****************************************************************************80
{
   int ndim = 5;     // dimension of space
   double beta = 1.; // parameter of Oscillatory function
   testfun tf(ndim,beta);

   tf.InternalIntegrator();

   return 0;
}
#endif

