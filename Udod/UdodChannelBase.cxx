#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <numeric>

// #include "UdodFinalState.h"

#include "UdodChannelBase.h"
#include "CkTable.h"
#include "Udod.h"
#include <cstdlib>

using namespace std;
using namespace UDOD;

//------------------------------------------------------------------------------
template <class Atype>
inline static Atype SQR(Atype x)
{
   return(x*x);
}

double UdodChannelBase::Random() {
//  return parent->Random();
  return 0;
}

//------------------------------------------------------------------------------
void UdodChannelBase::Initialize(std::vector<double>& Mass)
{
   npart = Mass.size();
   if( npart < 2 || npart > 20 ) {
      cout << "UdodChannel::UdodChannel ERROR"
           " Wrong number of decay particles: " << npart << endl;
      exit(1);
   }

   // ATTENTION: fist index of mass_v is 1 and mass_v[0] must be 0.
   mass_v.resize(npart+1,0.);
   for(unsigned int i = 0; i < npart; i++) {
      mass_v[i+1] = Mass[i];
   }

   // prepare vector for final state particles
//   parent->DFS.resize(npart);

   // Calculate Ck tables for KProcedure
   tbls.resize(npart,(CkTable*)0);
   for(unsigned int k = 2; k < npart; k++) {
      tbls[k] = new CkTable(k);
   }
}

//------------------------------------------------------------------------------
UdodChannelBase::~UdodChannelBase()
{
   // clean up CkTables
   vector<CkTable*>::iterator it = tbls.begin();
   for(; it != tbls.end(); it++) {
      if( *it ) {
         delete *it;
      }
   }
}

//------------------------------------------------------------------------------
void UdodChannelBase::CouplesProduction()
{
/*
   // Generate two particle decay X-> Particle^+ + Particle^-

   long double Ex  = parent->FS[XSTATE].Get_E();
   long double Px  = parent->FS[XSTATE].Get_P();
   long double Tx  = parent->FS[XSTATE].Get_Theta();
   long double Phx = parent->FS[XSTATE].Get_Phi();
   long double mx  = sqrtl(Ex * Ex - Px * Px);

   //----for particle with negative charge in X rest system
   long double EPart1 = 0.5 * mx;
   long double PPart1 = sqrtl(SQR(EPart1) - SQR(mass_v[1]));

   //  cout << "PPart1 = " << PPart1 << endl;

   long double ThetaPart1 = acosl(2. * Random() - 1.);
   long double PhiPart1   = 2.* M_PI * Random();

   //------------------------Lorenz Transformation--------------
   //------to Dekart system
   long double ip1 = PPart1 * sinl(ThetaPart1) * cosl(PhiPart1);
   long double jp1 = PPart1 * sinl(ThetaPart1) * sinl(PhiPart1);
   long double kp1 = PPart1 * cosl(ThetaPart1);
   //------to laboratorial system( ip & jp stay as they are)
   long double kp_lab1 = Ex / mx * kp1 + Px / mx * EPart1;
   long double E_lab1 =  Ex / mx * EPart1  + Px / mx * kp1;
   //------matrix for rotation
   long double a11=cosl(Tx)*cosl(Phx), a12=-sinl(Phx), a13=sinl(Tx)*cosl(Phx);
   long double a21=cosl(Tx)*sinl(Phx), a22=cosl(Phx),  a23=sinl(Tx)*sinl(Phx);
   long double a31=-sinl(Tx),          a32=0.,         a33=cosl(Tx);
   //------rotation
   long double px1 = a11 * ip1 + a12 * jp1 + a13 * kp_lab1;
   long double py1 = a21 * ip1 + a22 * jp1 + a23 * kp_lab1;
   long double pz1 = a31 * ip1 + a32 * jp1 + a33 * kp_lab1;
   //-----------------------------------------------------------

   //----for particle with positive charge in X rest system
   long double EPart2 = EPart1;
   long double PPart2, ThetaPart2, PhiPart2;
   //------------------------Lorenz Transformation--------------
   //------to Dekart system
   long double ip2 = -ip1;
   long double jp2 = -jp1;
   long double kp2 = -kp1;
   //------to laboratorial system
   long double kp_lab2 = Ex / mx * kp2 + Px / mx * EPart2;
   long double E_lab2 =  Ex / mx * EPart2  + Px / mx * kp2;
   //------rotation
   long double px2 = a11 * ip2 + a12 * jp2 + a13 * kp_lab2;
   long double py2 = a21 * ip2 + a22 * jp2 + a23 * kp_lab2;
   long double pz2 = a31 * ip2 + a32 * jp2 + a33 * kp_lab2;
   //-----------------------------------------------------------

   //--back to spherical system +
   //  correction of the angle phi of the all system
   long double PhiCorrector = 2. * M_PI * Random();

   //----for particle with negative charge
   EPart1     = E_lab1;
   PPart1     = sqrtl(px1 * px1 + py1 * py1 + pz1 * pz1);
   ThetaPart1 = atan2(sqrt(px1*px1 + py1*py1), pz1);
   PhiPart1   = atan2l(py1 , px1) + PhiCorrector;

   //----for particle with positive charge
   EPart2     = E_lab2;
   PPart2     = sqrtl(px2 * px2 + py2 * py2 + pz2 * pz2);
   ThetaPart2 = atan2(sqrt(px2*px2 + py2*py2), pz2);
   PhiPart2   = atan2l(py2 , px2) + PhiCorrector;

   //----setting
   parent->DFS[0].Set(EPart1,PPart1,ThetaPart1,PhiPart1);
   parent->DFS[1].Set(EPart2,PPart2,ThetaPart2,PhiPart2);

   parent->FS[ELECTRON].Set_Phi(parent->FS[ELECTRON].Get_Phi() + PhiCorrector);
   parent->FS[POSITRON].Set_Phi(parent->FS[POSITRON].Get_Phi() + PhiCorrector);
*/
}

//------------------------------------------------------------------------------
void UdodChannelBase::KProcedure()
{
/*   // Generate N-particles decay taking into account
   // conservation laws of energy and momentum.
   //
   // Notations and numbering made in accordance with the book:
   // Копылов Г.И. Основы кинематики резонансов, М., Наука, 1970 г.
   // (см. так же препринт ИЯФ 2000-78)
   // Note: instead of the tilde in variables, we use an underscore
   // on the end of the names.

   long double fs_array[3];

   //----k1: 4-momentum of N-particles in LAB system
   long double En = parent->FS[XSTATE].Get_E();
   long double Pn = parent->FS[XSTATE].Get_P();
   long double Pn_x = parent->FS[XSTATE].Get_Px();
   long double Pn_y = parent->FS[XSTATE].Get_Py();
   long double Pn_z = parent->FS[XSTATE].Get_Pz();
   long double Mn = En*En - Pn*Pn;
   if( Mn < 0 ) {
      cout << "UdodChannel::KProcedure ERROR"
           " Square of effective mass = " << Mn << endl;
      exit(1);
   }
   Mn = sqrtl(Mn);

   //----k2: Sum of the masses of particles (mu_n)
   long double mu_n = accumulate(mass_v.begin(), mass_v.end(),
                                 0., plus<double>() );
   if( mu_n < 0 ) {
      cout << "UdodChannel::KProcedure ERROR"
           " Total mass = " << mu_n << endl;
      exit(1);
   }

   //----k3: Energy release of N-particles (Tn_) in center of mass system
   long double Tn_ = Mn - mu_n;
   if( Tn_ < 0 ) {
      cout << "UdodChannel::KProcedure ERROR"
           " Energy release = " << Tn_ << endl;
      exit(1);
   }

   // Cycle for K = N, N-1, ..., 2 to generate momentum of
   // particle #K in the system of rest particles 1,...,K (X_k)
   // (denote them as the compound particle Ok)
   // and translate this momentum to LAB system.
   for(unsigned int k = npart; k > 1; k--) {
      //----k4: Generate energy release for k-1 particles (Tkm1_ = T_{k-1})
      long double Tkm1_ = 0;
      if( k > 2 ) {
         CkTable* tbl = tbls[k-1];
         long double ksi = tbl->GetKsi( Random() );
         if( ksi < 0. || ksi > 1. ) {
            cout << "UdodChannel::KProcedure ERROR"
                 " ksi = " << ksi << endl;
            exit(1);
         }
         Tkm1_ = Tn_ * ksi;
      }

      //----k5: Mass for k-1 particles
      mu_n -= mass_v[k];
      if( mu_n < 0 ) {
         cout << "UdodChannel::KProcedure ERROR"
              " mu_n(k=" << k << ")= " << mu_n << endl;
         exit(1);
      }

      //----k6: Effective mass of system of k-1 particles (Mkm1 = M_{k-1})
      long double Mkm1 = mu_n + Tkm1_;

      //----k7: Energy and momentum of #k particle in Ok rest system
      long double omega_k_ = (SQR(Mn) + SQR(mass_v[k]) - SQR(Mkm1)) / (2*Mn);
      if( omega_k_ < 0 ) {
         cout << "UdodChannel::KProcedure ERROR"
              "omega_{k=" << k << "} = " << omega_k_ << endl;
         exit(1);
      }
      long double Pk_ = SQR(omega_k_) - SQR(mass_v[k]);
      if( Pk_ < 0 ) {
         cout << "UdodChannel::KProcedure ERROR"
              "  SQR(P_{k=" << k << "}) = " << Pk_ << endl;
         exit (1);
      }
      Pk_ = sqrtl(Pk_);

      //----k8: Generate cosine theta for particle #k in the Ok rest system
      long double cos_theta = 2.*Random() - 1.;
      long double sin_theta = sqrtl(1.-SQR(cos_theta));

      //----k9: Generate phi angle for particle #k in the Ok rest system
      long double phi       = 2.*M_PI*Random();
      long double cos_phi   = cosl(phi);
      long double sin_phi   = sinl(phi);

      //----k10: Momentum of particle #k in the Ok rest system (Pk_x,y,z)
      long double Pk_x_ = Pk_*sin_theta*cos_phi;
      long double Pk_y_ = Pk_*sin_theta*sin_phi;
      long double Pk_z_ = Pk_*cos_theta;

      //----k11: Energy of particle #k in the LAB system
      long double sc_prod = Pn_x*Pk_x_+Pn_y*Pk_y_+Pn_z*Pk_z_;
      long double omega_k = (En*omega_k_ + sc_prod) / Mn;

      //----k12: Momentum of particle #k in the LAB system
      long double lambda = (omega_k + omega_k_) / (En + Mn);
      long double Pk_x = Pk_x_ + lambda * Pn_x;
      long double Pk_y = Pk_y_ + lambda * Pn_y;
      long double Pk_z = Pk_z_ + lambda * Pn_z;

      // Save energy and momentum of particle #k
      fs_array[0] = Pk_x;
      fs_array[1] = Pk_y;
      fs_array[2] = Pk_z;
      parent->DFS[k-1].Set(omega_k, fs_array);

      //----k13: 4-momentum of k-1 particles (new N=k-1)
      En   -= omega_k;
      Pn_x -= Pk_x;
      Pn_y -= Pk_y;
      Pn_z -= Pk_z;
      Pn    = sqrtl(Pn_x*Pn_x + Pn_y*Pn_y + Pn_z*Pn_z);
      Mn = Mkm1;

   } // end of for(unsigned int k...

   //----k14: Momentum and energy of the particle #1
   long double P1 = sqrtl(Pn_x*Pn_x + Pn_y*Pn_y + Pn_z*Pn_z);
   long double E1 = En;

   //----k15: Check that everything is correct
   long double M1 = sqrtl(E1*E1 - P1*E1);
   if( fabs(M1 - mass_v[1]) > 1.e-6 ) {
      cout << "UdodChannel::KProcedure ERROR"
           " M1 not equal to the mass of the particle #1" << endl
           << " M1 = " << M1 << " P1 = " << P1 << " E1 = " << E1 << endl;
      exit(1);
   }

   // Save energy and momentum of particle #1
   fs_array[0] = Pn_x;
   fs_array[1] = Pn_y;
   fs_array[2] = Pn_z;
   parent->DFS[0].Set(En, fs_array);
*/
   return;
}

void UdodChannelBase::ShowResults()
{
/*  cout << " Decay list: " << endl;
  for(unsigned int l = 0; l < npart; l++) {
    cout << " k# " << l
	 << " : E_k= " << setprecision(15) << parent->DFS[l].Get_E() 
	 << " P_k= "   << setprecision(15) << parent->DFS[l].Get_P() 
	 << " M_k= "   << setprecision(15) << parent->DFS[l].Get_M() << endl;
  }
*/
  cout << "it's prototype of function" << endl;
}
