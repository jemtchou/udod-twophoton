#ifndef UdodCut_h
#define UdodCut_h

#include <iostream>

namespace UDOD
{
   class UdodCut
   {
      public:
         UdodCut() {
            fMinVal = 0.;
            fMaxVal = 10e9;
            fName = "Unknown cut";
            fUnits = "";
         }

         UdodCut(long double min, long double max, 
                 const std::string name= "", const std::string units= "") {
            fMinVal = min;
            fMaxVal = max;
            fName = name;
            fUnits = units;
         }

         virtual ~UdodCut() {}

         // Access methods
         void SetCut(long double min, long double max, 
                     const std::string name="", const std::string units= "") {
            fMinVal = min;
            fMaxVal = max;
            fName = name;
            fUnits = units;
         }

         void SetMinimum(long double min) {
            fMinVal = min;
         }

         void SetMaximum(long double max) {
            fMaxVal = max;
         }

         long double Min() {
            return fMinVal;
         }

         long double Max() {
            return fMaxVal;
         }

         // Check variable against the cut
         bool Check(long double val) {
            return ( val>fMinVal && val<fMaxVal );
         }

         // Print cut bounds
         void Print() {
            std::cout << fName << "\t[" << fUnits << "] : min "
                      << fMinVal << "\t" << " max " << fMaxVal << std::endl;
         }

      private:
         long double fMinVal;
         long double fMaxVal;
         std::string fName;
         std::string fUnits;
   };
}

#endif
