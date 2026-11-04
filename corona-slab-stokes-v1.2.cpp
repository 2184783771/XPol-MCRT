#include <math.h>
#include "corona-slab-stokes-v1.2.h"
#include "corona-stokes-v1.2.h"
#include "specmath-stokes-v1.2.h"


void CoronaApplySlabDefaults(CoronaRegion& Corona,XRBSeedDisk& SeedDisk,
                             bool readCoronaSource,bool readSeedSource,
                             bool readSeedAngular,bool readDiskReflection)
{
  if (!readCoronaSource) Corona.source = false;
  if (!readSeedSource) SeedDisk.source = true;
  if (!readSeedAngular) SeedDisk.angular = 1;
  if (!readDiskReflection) SeedDisk.reflection = false;
}


void CoronaSampleSlab(CoronaRegion Corona,double r1,double r2,
                      double r3,double pos[])
{
  double R,phi;

  R = Corona.radius*sqrt(r1);
  phi = 2*pi*r2;

  pos[0] = Corona.pos[0]+R*cos(phi);
  pos[1] = Corona.pos[1]+R*sin(phi);
  pos[2] = Corona.pos[2]+Corona.height*(r3-0.5);
}


void CoronaFindSlabBorder(CoronaRegion Corona,double& tpre,
                          double pos[],double direct[])
{
  double h,r,t,LocalPos[3];

  t = tpre;
  r = Corona.radius;
  h = 0.5*Corona.height;
  VecDiff(LocalPos,pos,Corona.pos);

  if (r>0 && h>0)
  {
    CylinderCoinc(LocalPos,direct,t,r,h,0);
    RingCoinc(LocalPos,direct,t,0,r,-h);
    RingCoinc(LocalPos,direct,t,0,r,h);
  }

  if (t<tpre)
  {
    tpre = t;
  }
}


bool CoronaInsideSlab(CoronaRegion Corona,double pos[],double gap)
{
  double LocalPos[3],r,h,eps;

  VecDiff(LocalPos,pos,Corona.pos);
  r = Corona.radius;
  h = 0.5*Corona.height;
  eps = CoronaBoundaryEps(r,h,gap);

  return (r>0 && h>0 &&
          LocalPos[0]*LocalPos[0]+LocalPos[1]*LocalPos[1]<=(r+eps)*(r+eps) &&
          fabs(LocalPos[2])<=h+eps);
}
