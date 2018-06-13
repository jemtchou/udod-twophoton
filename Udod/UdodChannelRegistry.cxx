#include <cmath>
#include <iostream>
#include <vector>
#include <cassert>

#include "UdodChannelRegistry.h"
#include "Udod.h"
#include "UdodChannelBase.h"


using namespace std;
using namespace UDOD;

UdodChannelRegistry::UdodChannelRegistry(Udod* p)
{
  parent=p;
}

void UdodChannelRegistry::AddChannel(UdodChannelBase* newChannel)
{
   channels.push_back(newChannel);
}

UdodChannelRegistry::~UdodChannelRegistry()
{
   vector<UdodChannelBase*>::iterator it = channels.begin();
   for(; it != channels.end(); it++) {
      if( *it ) {
         delete *it;
      }
   }
}


void UdodChannelRegistry::Initialize()
{
  double currentSumXS = 0;
  for (unsigned int i=0; i<channels.size();i++)
    {
      if(channels[i]->IsActive())
	{
	  currentSumXS += channels[i]->GetIntegralCS();
	  Intervals.push_back(currentSumXS);
	  NInterval.push_back(i);
	}
    }
}

UdodChannelBase* UdodChannelRegistry::SelectChannel()
{
  assert(Intervals.size()>0);
  double a = parent->Random()*Intervals.back();
  for(unsigned int i=0; i<Intervals.size();i++)
    {
      if (a <= Intervals[i])
	{
	  return channels[NInterval[i]];
	}
    }
  return 0;
}
