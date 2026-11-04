#include "corona-lamppost-stokes-v1.2.h"
#include "corona-sphere-stokes-v1.2.h"


void CoronaApplyLampPostDefaults(CoronaRegion& Corona,XRBSeedDisk& SeedDisk,
                                 bool readCoronaSource,bool readSeedSource,
                                 bool readSeedAngular,bool readDiskReflection)
{
  if (!readCoronaSource) Corona.source = true;
  if (!readSeedSource) SeedDisk.source = false;
  if (!readDiskReflection) SeedDisk.reflection = true;
}


void CoronaSampleLampPost(CoronaRegion Corona,double r1,double r2,
                          double r3,double pos[])
{
  int i;

  if (Corona.radius>0)
  {
    CoronaSampleSphere(Corona,r1,r2,r3,pos);
  }
  else
  {
    for(i=0; i<=2; i++)
    {
      pos[i] = Corona.pos[i];
    }
  }
}


void CoronaFindLampPostBorder(CoronaRegion Corona,double& tpre,
                              double pos[],double direct[])
{
  if (Corona.radius>0)
  {
    CoronaFindSphereBorder(Corona,tpre,pos,direct);
  }
}


bool CoronaInsideLampPost(CoronaRegion Corona,double pos[],double gap)
{
  if (Corona.radius>0)
  {
    return CoronaInsideSphere(Corona,pos,gap);
  }
  return false;
}
