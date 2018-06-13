#include "ModelSigmaGGpi0.hh"

#define _USE_MATH_DEFINES //for M_PI
#include <cmath>  
#include <cstdlib> // for exit(1)
#include <complex>
#include <iostream>
#include "CkTable.h"
#include "Udod.h"
#include "UdodConstants.h"

using namespace std;
using namespace UDOD;


double ModelSigmaGGpi0::SigmaGG(double W2, double theta, double s, double s1, double s2, double t1, double t2)
{
  return alpha/W2*15.396*pow(10,6)/(4*M_PI*M_PI); //nBn
}

