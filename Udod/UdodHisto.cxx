#include "UdodHisto.h"

#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"

#include <iostream>
#include <cstdio>

using namespace UDOD;

UdodHisto* UdodHisto::fPtr = 0;

UdodHisto::UdodHisto()
{
  fout = new TFile("udod.root", "RECREATE");
  fout->cd();
  int i=0;
  for(i=0; i<NHIST; i++)
    {
      h1[i]=0;
      h2[i]=0;
    }
}

UdodHisto::~UdodHisto()
{
  int i=0;
   for(i=0; i<NHIST; i++)
    {
      if(h1[i]!=0) delete h1[i];
      if(h2[i]!=0) delete h2[i];
    }
   delete fout;
}

void UdodHisto::Hist1d(int id, std::string title, int nbinsx, double xmin, double xmax)
{
  if(id>=NHIST)
    {
      std::cerr << "Can book only " << NHIST << " histograms" << std::endl;
      return;
    }
  if(h1[id]!=0)
    {
      std::cerr << "Histogram " << id << " already booked" << std::endl;
      return;
    }
  char hname[20];
  sprintf(hname,"h1_%d",id);
  h1[id] = new TH1F(hname,title.c_str(),nbinsx,xmin,xmax);
}

void UdodHisto::Hist2d(int id, std::string title, int nbinsx, double xmin, double xmax, int nbinsy, double ymin, double ymax)
{
    if(id>=NHIST)
    {
      std::cerr << "Can book only " << NHIST << " histograms" << std::endl;
      return;
    }
  if(h2[id]!=0)
    {
      std::cerr << "Histogram " << id << " already booked" << std::endl;
      return;
    }
  char hname[20];
  sprintf(hname,"h2_%d",id);
  h2[id] = new TH2F(hname,title.c_str(),nbinsx,xmin,xmax,nbinsy,ymin,ymax);
}

void UdodHisto::Fill(int id, double x, double w)
{
  if(id>=NHIST)
    {
      std::cerr << "Histogram " << id << " doesn't exist" << std::endl;
      return;
    }
  if(h1[id]==0)
    {
      std::cerr << "Histogram " << id << " is not booked" << std::endl;
      return;
    }
  h1[id]->Fill(x,w);
}

void UdodHisto::Fill(int id, double x, double y, double w)
{
  if(id>=NHIST)
    {
      std::cerr << "Histogram " << id << " doesn't exist" << std::endl;
      return;
    }
  if(h2[id]==0)
    {
      std::cerr << "Histogram " << id << " is not booked" << std::endl;
      return;
    }
  h2[id]->Fill(x,y,w);
}

void UdodHisto::Save()
{
  fout->cd();
  int i=0;
  for(i=0; i<NHIST; i++)
    {
      if(h1[i]!=0) h1[i]->Write();
      if(h2[i]!=0) h2[i]->Write();
    }
  fout->Print();
  fout->Close();
}
