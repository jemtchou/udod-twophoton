#ifndef VegasRandom_h
#define VegasRandom_h

class VegasRandom
{
 public:
  VegasRandom(int ndim) { ndim_ = ndim; };

  // IniRandom sets up the random-number generator to produce a
  // sequence of at least n ndim_-dimensional random vectors.
  virtual void IniRandom(const long int n) = 0;
  
  // GetRandom retrieves one random vector.
  virtual void GetRandom(double *x) = 0;
  
  // SkipRandom skips over n random vectors.
  virtual void SkipRandom(long int n) = 0;

 protected:
  int ndim_;

};

#endif
