#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <numeric>

#include "rho_pole.hh"

#include "CkTable.h"
#include "Udod.h"
#include "UdodConstants.h"

using namespace std;
using namespace UDOD;
 
double rho_pole::hT(double Q2)
  {
    return (mRhoSq/(mRhoSq+Q2))*(mRhoSq/(mRhoSq+Q2));
  }

  double rho_pole::hS(double Q2)
  {
    return ksi*Q2/mRhoSq*(mRhoSq/(mRhoSq+Q2))*(mRhoSq/(mRhoSq+Q2));
  }
 
