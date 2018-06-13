#ifndef VegasUtils_h
#define VegasUtils_h

class VegasUtils
{
 public:
  static double Sq(const double x);

  static double Min(const double a, const double b);

  static double Max(const double a, const double b);

  static double Weight(const double sum, const double sqsum, const long int n);

  // (a < 0) ? -1 : 0 
  static int NegQ(int a);

  // (a < 0) ? 0 : a 
  static int IDim(int a);

  // (a < b) ? a : b 
  static int IMin(int a, int b);

  // (a > b) ? a : b 
  static int IMax(int a, int b);

  //(a == 0) ? 0 : -1 
  static int TrueQ(int a);

  // a + (a == 0) 
  static int Min1(int a);

  // abs(a) + (a == 0) 
  static int Abs1(int a);
  
  //W.J. Kennedy and J.E. Gentle, Statistical computing, p. 116  
  static double Normal(const double x);
  
  static double ChiSquare(const double x, const int df);

  //	Gaussian error function
  //	= 2/Sqrt[Pi] Integrate[Exp[-t^2], {t, 0, x}]
  //	Code from Takuya Ooura's gamerf2a.f
  //	http://www.kurims.kyoto-u.ac.jp/~ooura/gamerf.html
  static double Erfc(const double x);

  static double Erf(const double x);

 private:
  static const double NOTZERO=0x1p-104;

};

#endif
