#include <iostream>
#include <cmath>
#include <algorithm>
#include <functional>

#include "CkTable.h"

using namespace std;

//------------------------------------------------------------------------------
CkTable::CkTable(int npart, int nsteps)
{
   // check values:
   if( npart < 3 ) {
      cout << "ERROR::CkTable::CkTable() number of particles npart="
           << npart << " must be greater than 3" << endl;
      exit(1);
   }
   if( nsteps < 20 ) {
      cout << "ERROR::CkTable::CkTable() number of steps nsteps="
           << nsteps << " must be greater than 20" << endl;
      exit(1);
   }

   // fill tables
   nparticles    = npart;
   full_integral = CkCalculation(nsteps);
}

//------------------------------------------------------------------------------
double CkTable::GetKsi(double Ck)
{
   // binary search in the table
   vector<double>::iterator p =
      lower_bound(ck_table.begin(), ck_table.end(), Ck);

   if( p == ck_table.end() ) {
      cout << "WARNING::CkTable :: GetKsi: "
           << Ck << " not found in ck_table" << endl;
      return 1.;
   }

   int idx = distance(ck_table.begin(), p);
   // return LinInterpol(Ck, idx);
   return QuadInterpol(Ck, idx);
}

//------------------------------------------------------------------------------
double CkTable::Integrand(double ksi)
{
   // calculate function x^(3k-5)/2 * sqrt(1-x)
   static const double eps = 2.2e-16;

   double x = fabs(ksi);
   double y = fabs(1.-ksi);
   if( x*y < eps ) {
      return 0;
   }
   div_t k_div = div(3*nparticles-5, 2);
   double f = pow_n(x,k_div.quot);
   if( k_div.rem ) {
      f *= sqrt(x*y);
   } else {
      f *= sqrt(y);
   }

   return f;
}

//------------------------------------------------------------------------------
double CkTable::CkCalculation(int nsteps)
{
   step = 1./double(nsteps);
   unsigned int tbl_size = nsteps+1;
   ksi_table.resize(tbl_size);
   ck_table.resize(tbl_size);

   // calculate the integrand function on the mesh
   // and temporary save results in ck_table
   double ksi = 0.;
   int i      = 0;
   for(; i < nsteps; i++) {
      ksi_table[i] = ksi;
      ck_table[i]  = Integrand(ksi);
      ksi          += step;
   }
   ksi          = 1.; // the last point is exactly 1
   ksi_table[i] = ksi;
   ck_table[i]  = Integrand(ksi);

   // Use Simpson formula for calculating integrals in a bounds [0,ksi_i]
   // (actually these are multiplied to 6/step)
   vector<double> cumulative(tbl_size,0.);
   for(unsigned int i = 1; i < tbl_size; i++) {
      ksi = 0.5*(ksi_table[i-1] + ksi_table[i]);
      cumulative[i] = cumulative[i-1] +
                      (ck_table[i-1] + 4*Integrand(ksi) + ck_table[i]);
   }

   // normalize cumulative sums such that cumulative[nsteps] == 1
   // and copy them in ck_table
   double norm = 1./cumulative[nsteps];
   transform(cumulative.begin(), cumulative.end(), // from
             ck_table.begin(),                     // to
             bind2nd(multiplies<double>(),norm) ); // action

   // calculate full_integral and return it
   double ret = cumulative[nsteps] * step/6;
   return ret;
}

//------------------------------------------------------------------------------
double CkTable::LinInterpol(double ck, unsigned int up_index)
{
   // reverse liner interpolation in table
   // (find ksi corresponding to ck)
   if(up_index == 0) {
      return 0.;
   }

   double ck0 = ck_table[up_index - 1];
   double ck1 = ck_table[up_index];
   double p = (ck-ck0)/(ck1-ck0);

   double ksi0 = ksi_table[up_index - 1];
   return ksi0 + p*step;
}

//------------------------------------------------------------------------------
double CkTable::QuadInterpol(double ck, unsigned int up_index)
{
   // reverse quadratic interpolation in table
   // (find ksi corresponding to ck)
   if(up_index == 0) {
      return 0.;
   }

   double ckM = ck_table[up_index - 1];
   double ckO = ck_table[up_index];
   double ckP = 0;
   double ksiO = ksi_table[up_index];
   if( up_index+1 < ck_table.size() ) {
      ckP = ck_table[up_index + 1];
   } else {
      // this is case of the last sell in table
      // cout << " last sell in table:" << endl;
      ckP = ckO;
      ckO = ckM;
      ckM = ck_table[up_index - 2];
      ksiO = ksi_table[up_index - 1];
   }

   double A = ckP - 2*ckO + ckM;
   double B = ckP - ckM;
   double C = 2*(ckO - ck);

   double D = B*B - 4*A*C;
   if( D < 0 ) {
      cout << " ERROR: CkTable::QuadInterpol D= " << D << endl;
      return LinInterpol(ck,up_index);
   }
   D=sqrt(D);

   // we know that ckP > ckO > ckM => B > 0 => we should use (-B+D) solution
   double p = (-B + D)/(2*A);
   if( fabs(p) > 1 ) {
      cout << " ERROR: CkTable::QuadInterpol p(+)= " << p
           << " p(-)= " << (-B - D)/(2*A) << endl;
      return LinInterpol(ck,up_index);
   }

   return ksiO + p*step;
}

//------------------------------------------------------------------------------
double CkTable::pow_n(double x, unsigned int n)
{
   double res = 1.;
   if( n != 0 ) {
      for(unsigned int i = 0x1; ; i <<= 1) { /* cycle for bits */
         if( n & i ) {                       /* check bit-i in n */
            res *= x;
            if( (n &= ~i) == 0 ) {           /* set bit-i to 0 */
               break;
            }
         }
         x *= x;
      }
   }
   return res;
}

//------------------------------------------------------------------------------
#ifdef MAINCKTABLE
#include <cstdio>

//------------------------------------------------------------------------------
double AnalyticalInt(unsigned int n)
{
   // Integral calculated analytically:
   //   \int_0^1 dx [x^{(3n-5)/2} \sqrt{1-x}] =
   //            \sqrt(\pi)/2 * Γ((3n-3)/2)/Γ(3n/2)
   //   n=2 -> pi/8            n=3 -> 16/105
   //   n=4 -> 7pi/256         n=5 -> 512/9009
   //                       (n-1)(3n-1)(3n+1)
   //   Int(n+2) = Int(n) x -----------------;
   //                         n  (3n+2)(3n+4)

   static vector<double> Save(20,1.);
   if( n < 2 || n >= Save.size() ) {
      return 0;
   }

   if( Save[0] > 0 ) {
      Save[0] = Save[1] = 0.;
      Save[2] = M_PI/8;
      Save[3] = 16./105;
      for(unsigned int k = 2; k < Save.size()-2; k++) {
         Save[k+2] = Save[k] * (k-1)*(3*k-1)*(3*k+1)/double(k*(3*k+2)*(3*k+4));
      }
   }
   return Save[n];
}

//------------------------------------------------------------------------------
double AnalyticalFun3(double x)
{
   // analytically calculated function Ck(k=3)
   double y = 1 - x;
   return 1.-0.125*y*sqrt(y)*(35.+y*(-42.+y*15.));
}

//------------------------------------------------------------------------------
double AnalyticalFun4(double x)
{
   // analytically calculated function Ck(k=4)
   double y = 1 - x;
   return M_2_PI * ( asin(sqrt(x)) +
                     sqrt(x*y)/105.*(-105.+x*(-70.+x*(-56.+x*(-48.+x*384.)))) );
}

//------------------------------------------------------------------------------
double AnalyticalFun5(double x)
{
   // analytically calculated function Ck(k=5)
   double y = 1 - x;
   return 1.-3003./256*y*sqrt(y)*
          (1.+y*(-3.+y*(30./7+y*(-10./3+y*(15./11+y*(-3./13))))));
}

//------------------------------------------------------------------------------
double AnalyticalFun(int k, double x)
{
   double ret = 0;
   switch(k) {
   case 3: ret = AnalyticalFun3(x); break;
   case 4: ret = AnalyticalFun4(x); break;
   case 5: ret = AnalyticalFun5(x); break;
   }
   return ret;
}

//------------------------------------------------------------------------------
int main()
{
   printf("\n +++ Test of internal integral +++\n");
   printf("   Npar    TestInt                 AnalyticalInt      (eps)\n");
   for(int npar = 3; npar < 20; npar++) {
      CkTable test(npar);
      double Ck_int = test.GetFullIntegral();
      double A_int = AnalyticalInt(npar);
      double eps = fabs(Ck_int/A_int - 1.);
      printf(" %3i %23.15f %23.15f  (%8.2e)\n",
             test.GetNparticles(), Ck_int, A_int, eps);
   }

   printf("\n\n +++ Test search in table function +++\n");
   int npar = 5;
   CkTable test(npar);
   double eps_max = 0.;
   printf(" ck        ksi        a_ck      (eps)\n");
   for(double ck = 0.001; ck < 1.; ck += 0.005) {
      double ksi = test.GetKsi(ck);
      double a_ck = AnalyticalFun(npar,ksi);
      double eps = fabs(a_ck/ck - 1.);
      if( eps > eps_max ) {
         eps_max = eps;
      }
      printf(" %5.3f %10.4f %12.6f  (%8.2e)\n", ck, ksi, a_ck, eps);
   }
   printf(" -------------------------------\n"
          " Max eps(npar=%i)= %8.2e \n", npar, eps_max);

   return 0;
}
#endif

