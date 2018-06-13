#include "VegasUtils.h"
#include <cmath>

double VegasUtils::Sq(const double x){ return x*x; }

double VegasUtils::Min(const double a, const double b){ return (a < b) ? a : b;}

double VegasUtils::Max(const double a, const double b){ return (a > b) ? a : b;}

double VegasUtils::Weight(const double sum, const double sqsum, const long int n)
{
  const double w = sqrt(sqsum*n);
  return (n - 1)/Max((w + sum)*(w - sum), NOTZERO);
}

int VegasUtils::NegQ(int a) { return ((a) >> (sizeof(a)*8 - 1)); }

int VegasUtils::IDim(int a) { return ((a) & NegQ(-(a))); }

int VegasUtils::IMin(int a, int b) { return ((a) - IDim((a) - (b))); }

int VegasUtils::IMax(int a, int b) { return ((b) + IDim((a) - (b))); }

int VegasUtils::TrueQ(int a) { return NegQ((a) | (-a)); }

int VegasUtils::Min1(int a) { return ((a) + 1 + TrueQ(a)); }

int VegasUtils::Abs1(int a) { return (((a) ^ NegQ(a)) - NegQ((a) - 1)); }

double VegasUtils::Normal(const double x)
{
  return .5*erf(x/1.414213562373095048801689) + .5;
}

double VegasUtils::ChiSquare(const double x, const int df)
{
  double y;

  if( df <= 0 ) return -999;

  if( x <= 0 ) return 0;
  if( x > 1000*df ) return 1;

  if( df > 1000 ) {
    if( x < 2 ) return 0;
    y = 2./(9*df);
    y = (pow(x/df, 1/3.) - (1 - y))/sqrt(y);
    if( y > 5 ) return 1;
    if( y < -18.8055 ) return 0;
    return Normal(y);
  }

  y = .5*x;

  if( df & 1 ) {
    const double sqrty = sqrt(y);
    double h = erf(sqrty);
    int i;

    if( df == 1 ) return h;

    y = sqrty*exp(-y)/.8862269254527579825931;
    for( i = 3; i < df; i += 2 ) {
      h -= y;
      y *= x/i;
    }
    y = h - y;
  }
  else {
    double term = exp(-y), sum = term;
    int i;

    for( i = 1; i < df/2; ++i )
      sum += term *= y/i;
    y = 1 - sum;
  }

  return Max(0., y);
}

double Erfc(const double x)
{
  static const double c[] = {
    2.96316885199227378e-01, 6.12158644495538758e-02,
    1.81581125134637070e-01, 5.50942780056002085e-01,
    6.81866451424939493e-02, 1.53039662058770397e+00,
    1.56907543161966709e-02, 2.99957952311300634e+00,
    2.21290116681517573e-03, 4.95867777128246701e+00,
    1.91395813098742864e-04, 7.41471251099335407e+00,
    9.71013284010551623e-06, 1.04765104356545238e+01,
    1.66642447174307753e-07, 1.48455557345597957e+01,
    6.10399733098688199e+00, 1.26974899965115684e+01 };
  double y = x*x;
  y = exp(-y)*x*(
    c[0]/(y + c[1]) + c[2]/(y + c[3]) +
    c[4]/(y + c[5]) + c[6]/(y + c[7]) +
    c[8]/(y + c[9]) + c[10]/(y + c[11]) +
    c[12]/(y + c[13]) + c[14]/(y + c[15]) );
  if( x < c[16] ) y += 2/(exp(c[17]*x) + 1);
  return y;
}

double Erf(const double x)
{
  static const double c[] = {
    1.12837916709551257e+00,
    -3.76126389031833602e-01,
    1.12837916706621301e-01,
    -2.68661698447642378e-02,
    5.22387877685618101e-03,
    -8.49202435186918470e-04 };
  double y = fabs(x);
  if( y > .125 ) {
    y = 1 - Erfc(y);
    return (x > 0) ? y : -y;
  }
  y *= y;
  return x*(c[0] + y*(c[1] + y*(c[2] +
				y*(c[3] + y*(c[4] + y*c[5])))));
}

