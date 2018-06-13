#ifndef UdodChannelRegistry_h
#define UdodChannelRegistry_h

#include <iostream>
#include <vector>
#include "UdodChannelBase.h"

namespace UDOD{

  class Udod;

class UdodChannelRegistry
{
 public:
  UdodChannelRegistry(Udod* p);
  ~UdodChannelRegistry();

  void Initialize();

  void AddChannel(UdodChannelBase* newChannel);
  
  unsigned int GetNChannels() {return channels.size();}
  UdodChannelBase* GetChannel(int i) { return channels[i];}
  
  UdodChannelBase* SelectChannel();

 private:
  std::vector<UdodChannelBase*> channels;

  std::vector<double> Intervals;
  std::vector<int>    NInterval;

  Udod* parent;
};
}

#endif
