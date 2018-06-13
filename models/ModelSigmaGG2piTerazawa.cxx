#include "ModelSigmaGG2piTerazawa.hh"

#define _USE_MATH_DEFINES //for M_PI
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <numeric>

#include "CkTable.h"
#include "Udod.h"
#include "UdodConstants.h"

using namespace std;
using namespace UDOD;

double ModelSigmaGG2piTerazawa::fpi(double x)
{
  return pow(mV,4)/((x+mV*mV)*(1./2*x+mV*mV));			    
}
double ModelSigmaGG2piTerazawa::SigmaGG(double W2, double theta, double s, double s1, double s2, double t1, double t2)
{
  return M_PI*alpha*alpha/W2*pow((1-4*mpi*mpi/W2),1./2)*fpi(1./2*W2)*fpi(1./2*W2);
  //*15.396*pow(10,6)/(4*M_PI*M_PI);//nBn			    
}
