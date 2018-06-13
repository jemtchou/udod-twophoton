#define _USE_MATH_DEFINES //for M_PI

#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <numeric> 
#include <cstdlib> // for exit(1) 
#include "VMDc.hh"
#include "CkTable.h"
#include "Udod.h"
#include "UdodConstants.h"

using namespace std;
using namespace UDOD;

double VMDc::hT(double Q2)
{
    return rRho*(mRhoSq/(mRhoSq+Q2))*(mRhoSq/(mRhoSq+Q2))+
      rOmega*(mOmegaSq/(mOmegaSq+Q2))*(mOmegaSq/(mOmegaSq+Q2))+
      rPhi*(mPhiSq/(mPhiSq+Q2))*(mPhiSq/(mPhiSq+Q2))+
      rc*(m0sq/(m0sq+Q2));
}

double VMDc::hS(double Q2)
{
  return ksi*Q2*( 1/mRhoSq*rRho*(mRhoSq/(mRhoSq+Q2))*(mRhoSq/(mRhoSq+Q2))+
		  1/mOmegaSq*rOmega*(mOmegaSq/(mOmegaSq+Q2))*(mOmegaSq/(mOmegaSq+Q2))+
		  1/mPhiSq*rPhi*(mPhiSq/(mPhiSq+Q2))*(mPhiSq/(mPhiSq+Q2)) );
}
