#include <math.h>
#include "corona-stokes-v1.2.h"


void CoronaApplySourceDefaults(CoronaRegion& Corona,XRBSeedDisk& SeedDisk,
                               bool readCoronaSource,bool readSeedSource,
                               bool readSeedAngular,bool readDiskReflection)
{
  switch(Corona.type)
  {
    case 0:
      CoronaApplySlabDefaults(Corona,SeedDisk,readCoronaSource,
                              readSeedSource,readSeedAngular,
                              readDiskReflection);
      break;

    case 1:
      CoronaApplyToroidalDefaults(Corona,SeedDisk,readCoronaSource,
                                  readSeedSource,readSeedAngular,
                                  readDiskReflection);
      break;

    default:
      CoronaApplyLampPostDefaults(Corona,SeedDisk,readCoronaSource,
                                  readSeedSource,readSeedAngular,
                                  readDiskReflection);
      break;
  }
}


void CoronaSamplePosition(CoronaRegion Corona,double r1,double r2,
                          double r3,double pos[])
{
  switch(Corona.type)
  {
    case 0:
      CoronaSampleSlab(Corona,r1,r2,r3,pos);
      break;

    case 1:
      CoronaSampleToroidal(Corona,r1,r2,r3,pos);
      break;

    default:
      CoronaSampleLampPost(Corona,r1,r2,r3,pos);
      break;
  }
}


// sample a photon energy in keV from the corona power-law spectrum
double CoronaSampleEnergy(CoronaRegion Corona,double r)
{
  double emin,emax,gamma,min,max,energy;

  emin = Corona.emin;
  emax = Corona.emax;
  gamma = Corona.gamma;

  if (emin<=0 || emax<=0)
  {
    return 0;
  }
  if (emax<emin)
  {
    energy = emin;
    emin = emax;
    emax = energy;
  }
  if (emax==emin)
  {
    return emin;
  }

  if (gamma==1)
  {
    energy = emin*pow(emax/emin,r);
  }
  else
  {
    min = pow(emin,1-gamma);
    max = pow(emax,1-gamma);
    energy = pow(min+r*(max-min),1/(1-gamma));
  }
  return energy;
}


// if the corona is in the photon's flight direction, set tpre to the
// distance to the nearest corona border, otherwise leave it unchanged
void CoronaFindNextBorder(CoronaRegion Corona,double& tpre,
                          double pos[],double direct[])
{
  switch(Corona.type)
  {
    case 0:
      CoronaFindSlabBorder(Corona,tpre,pos,direct);
      break;

    case 1:
      CoronaFindToroidalBorder(Corona,tpre,pos,direct);
      break;

    default:
      CoronaFindLampPostBorder(Corona,tpre,pos,direct);
      break;
  }
}


bool CoronaInside(CoronaRegion Corona,double pos[],double gap)
{
  switch(Corona.type)
  {
    case 0:
      return CoronaInsideSlab(Corona,pos,gap);

    case 1:
      return CoronaInsideToroidal(Corona,pos,gap);

    default:
      return CoronaInsideLampPost(Corona,pos,gap);
  }
}


double CoronaBoundaryEps(double radius,double height,double gap)
{
  double scale,eps;

  scale = radius;
  if (height>scale)
  {
    scale = height;
  }
  if (scale<=0)
  {
    scale = 1;
  }

  if (gap>0)
  {
    eps = gap;
    if (eps<1.0e-12*scale)
    {
      eps = 1.0e-12*scale;
    }
  }
  else
  {
    eps = 0;
  }

  return eps;
}
