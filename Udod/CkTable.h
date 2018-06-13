#ifndef CKTABLE_H
#define CKTABLE_H

#include <vector>

class CkTable
{
   public:
      // default ctor: set number of particles and number of steps
      CkTable(int npart = 3, int nsteps = 200);

      // access functions:
      int GetNparticles()      {
         return nparticles;
      }
      double GetStep()         {
         return step;
      }
      double GetFullIntegral() {
         return full_integral;
      }

      // find ksi corresponding to given ck:
      double GetKsi(double Ck);

   private:
      int    nparticles;    // number of particles
      double step;          // number of interval of integration
      double full_integral; // integral on [0,1]

      std::vector<double> ksi_table; // ksi mesh
      std::vector<double> ck_table;  // Ck on this mesh

      // auxiliary functions:

      // this function integrated:
      double Integrand(double ksi);

      // major function to calculate Ck table:
      double CkCalculation(int nsteps);

      // interpolation in Ck table: find ksi corresponding to ck
      double LinInterpol(double ck, unsigned int up_index);
      double QuadInterpol(double ck, unsigned int up_index);

      // raise x to the positive integer power of n
      double pow_n(double x, unsigned int n);
};

#endif
