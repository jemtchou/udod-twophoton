#define _USE_MATH_DEFINES //for M_PI

#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <numeric>  
#include "rho_pole0.hh"
#include "CkTable.h"
#include "Udod.h"
#include "UdodConstants.h"

using namespace std;
using namespace UDOD;


double rho_pole0::hT(double Q2)
{
  return (mRhoSq/(mRhoSq+Q2))*(mRhoSq/(mRhoSq+Q2));
}

double rho_pole0::hS(double Q2)
{
  return 0;
}
  
