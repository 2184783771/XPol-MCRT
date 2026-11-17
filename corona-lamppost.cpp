#include "corona-lamppost.h"
#include "corona.h"
#include "specmath.h"
#include <math.h>

static void CoronaSampleSphere(CoronaRegion Corona,double r1,double r2,
                        double r3,double pos[])
{
  int i;
  double R,phi,theta,unit[3];

  R = Corona.radius*pow(r1,1.0/3);
  phi = 2*pi*r2;
  theta = acos(1-2*r3);
  CartCoord(unit,phi,theta);

  for(i=0; i<=2; i++)
  {
    pos[i] = Corona.pos[i]+R*unit[i];
  }
}


static void CoronaFindSphereBorder(CoronaRegion Corona,double& tpre,
                            double pos[],double direct[])
{
  double r,t;

  t = tpre;
  r = Corona.radius;

  if (r>0)
  {
    SphereSegCoinc(pos,direct,Corona.pos,t,r,pi/2,true);
  }

  if (t<tpre)
  {
    tpre = t;
  }
}


static bool CoronaInsideSphere(CoronaRegion Corona,double pos[],double gap)
{
  double LocalPos[3],r,eps;

  VecDiff(LocalPos,pos,Corona.pos);
  r = Corona.radius;
  eps = CoronaBoundaryEps(r,0,gap);

  return (r>0 && SProd(LocalPos,LocalPos)<=(r+eps)*(r+eps));
}



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
