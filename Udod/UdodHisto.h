#ifndef UdodHisto_h
#define UdodHisto_h

#include <string>
#include <iostream>

class TH2F;
class TH1F;
class TFile;

namespace UDOD{

#define NHIST 1000

class UdodHisto
{
 public:
  static UdodHisto* GetInstance()
    {
      if(! fPtr) fPtr = new UdodHisto();
      return fPtr;
    };
  ~UdodHisto();

  void Hist1d(int id, std::string title, int nbinsx, double xmin, double xmax);
  void Hist2d(int id, std::string title, int nbinsx, double xmin, double xmax, int nbinsy, double ymin, double ymax);
  void Fill(int id, double x, double w);
  void Fill(int id, double x, double y, double w);

  void Save();

 private:
  UdodHisto();

  TH1F* h1[NHIST];
  TH2F* h2[NHIST];

  TFile* fout;

  static UdodHisto* fPtr;
};
}

#endif
