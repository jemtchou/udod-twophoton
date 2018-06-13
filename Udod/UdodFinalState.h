#ifndef UdodFinalState_h
#define UdodFinalState_h

#include <cmath>
//#include "HepMC/SimpleVector.h"

namespace UDOD
{
   class UdodFinalState
   {
      public:
         UdodFinalState() {
            fE=0.;
            fP=0.;
            fTheta=0.;
            fPhi=0.;
         }
         ~UdodFinalState() {};

         void Set(long double e, long double p, 
                  long double theta, long double phi) {
            fE=e;
            fP=p;
            fTheta=theta;
            fPhi = NormalizePhi(phi);
         }

         void Set(long double e, const long double Pxyz[]) {
            fE=e;
            Set_Pvec(Pxyz);
         }

         long double E() const         {
            return fE;
         }
         long double P() const         {
            return fP;
         }
         long double Theta() const     {
            return fTheta;
         }
         long double Phi() const       {
            return fPhi;
         }
         long double Get_E() const     {
            return fE;
         }
         long double Get_P() const     {
            return fP;
         }
	 long double Get_M() const     {
	    return sqrt(fE*fE-fP*fP);
	 }
         long double Get_Theta() const {
            return fTheta;
         }
         long double Get_Phi() const   {
            return fPhi;
         }

         long double Get_Px() const    {
            return fP*sinl(fTheta)*cosl(fPhi);
         }
         long double Get_Py() const    {
            return fP*sinl(fTheta)*sinl(fPhi);
         }
         long double Get_Pz() const    {
            return fP*cosl(fTheta);
         }

         void Get_Pvec(long double Pxyz[]) const {
            long double Pt = fP*sinl(fTheta);
            Pxyz[0] = Pt*cosl(fPhi);        // Px
            Pxyz[1] = Pt*sinl(fPhi);        // Py
            Pxyz[2] = fP*cosl(fTheta);      // Pz
         }

         void Set_E(long double evar)     {
            fE = evar;
         }
         void Set_P(long double pvar)     {
            fP = pvar;
         }
         void Set_Theta(long double tvar) {
            fTheta = tvar;
         }
         void Set_Phi(long double phvar)  {
            fPhi = NormalizePhi(phvar);
         }
         void Set_Pvec(const long double Pxyz[]) {
            long double Px = Pxyz[0];
            long double Py = Pxyz[1];
            long double Pz = Pxyz[2];
            fP      = sqrtl( Px*Px + Py*Py + Pz*Pz );
            fTheta  = acosl( Px/fP );
            fPhi    = atan2l( Py,Px ); // range [-pi, pi]
         }

         long double NormalizePhi(long double angle) {
            while(angle > M_PI) {
               angle -= 2*M_PI;
            }
            while(angle < -M_PI) {
               angle += 2*M_PI;
            }
            return angle;
         }

      //   HepMC::FourVector GetFourVector() {
      //      return HepMC::FourVector(Get_Px(),Get_Py(),Get_Pz(),fE);
      //   }

      private:
         long double fE;
         long double fP;
         long double fTheta;
         long double fPhi;
   };
};

#endif
