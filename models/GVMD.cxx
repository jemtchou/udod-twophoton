#define _USE_MATH_DEFINES //for M_PI

#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <numeric> 
#include <cstdlib> // for exit(1)

#include "GVMD.hh"
#include "CkTable.h"
#include "Udod.h"
#include "UdodConstants.h"

using namespace std;
using namespace UDOD;
/*
double GVMD::P1(double Q2)
{
  return 1+Q2/m1sq;
}

double GVMD::P2(double Q2)
{
  return 1+Q2/m2sq;
}

double GVMD::hT(double Q2)
{
  return r*pow(P1(Q2),-2)+(1-r)*pow(P2(Q2),-1);
}

double GVMD::hS(double Q2)
{
   return ksi*(r*Q2/m1sq*pow(P1(Q2),-2)+(1-r)*(m2sq/Q2*log(P2(Q2))-pow(P2(Q2),-1)));
}
*/

double GVMD::hT(double Q2)
{
  return r*pow(P1,-2)+(1-r)*pow(P2,-1);
}

double GVMD::hS(double Q2)
{
   return ksi*(r*Q2/m1sq*pow(P1,-2)+(1-r)*(m2sq/Q2*log(P2)-pow(P2,-1)));
}
