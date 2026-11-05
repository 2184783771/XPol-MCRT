#include <math.h>
#include "corona-toroidal-stokes-v1.2.h"
#include "corona-stokes-v1.2.h"
#include "specmath-stokes-v1.2.h"

static double CoronaToroidalTubeRadius(CoronaRegion Corona)
{
  double tube;

  tube = Corona.height;
  if (tube<=0 && Corona.radius>0)
  {
    tube = 0.5*Corona.radius;
  }
  if (Corona.radius>0 && tube>=Corona.radius)
  {
    tube = 0.95*Corona.radius;
  }

  return tube;
}


static double SampleToroidalAngle(double major,double minor,double r)
{
  int i;
  double theta,k,f,df;

  theta = 2*pi*r;
  if (major<=0 || minor<=0)
  {
    return theta;
  }

  k = 2*minor/(3*major);
  for(i=0; i<8; i++)
  {
    f = theta+k*sin(theta)-2*pi*r;
    df = 1+k*cos(theta);
    if (df==0)
    {
      break;
    }
    theta -= f/df;
  }
  if (theta<0) theta = 0;
  if (theta>2*pi) theta = 2*pi;

  return theta;
}


static double SampleToroidalMinorRadius(double major,double minor,
                                        double theta,double r)
{
  int i;
  double low,high,mid,target,value,costheta;

  costheta = cos(theta);
  target = r*(0.5*major*minor*minor+
              costheta*minor*minor*minor/3.0);
  low = 0;
  high = minor;
  for(i=0; i<32; i++)
  {
    mid = 0.5*(low+high);
    value = 0.5*major*mid*mid+costheta*mid*mid*mid/3.0;
    if (value<target)
    {
      low = mid;
    }
    else
    {
      high = mid;
    }
  }

  return 0.5*(low+high);
}


void CoronaApplyToroidalDefaults(CoronaRegion& Corona,XRBSeedDisk& SeedDisk,
                                 bool readCoronaSource,bool readSeedSource,
                                 bool readSeedAngular,bool readDiskReflection)
{
  if (!readCoronaSource) Corona.source = false;
  if (!readSeedSource) SeedDisk.source = true;
  if (!readDiskReflection) SeedDisk.reflection = true;
  if (Corona.height<=0 && Corona.radius>0)
  {
    Corona.height = 0.5*Corona.radius;
  }
}


void CoronaSampleToroidal(CoronaRegion Corona,double r1,double r2,
                          double r3,double pos[])
{
  double major,minor,phi,theta,r,rho;

  major = Corona.radius;
  minor = CoronaToroidalTubeRadius(Corona);
  phi = 2*pi*r1;
  theta = SampleToroidalAngle(major,minor,r2);
  r = SampleToroidalMinorRadius(major,minor,theta,r3);
  rho = major+r*cos(theta);

  pos[0] = Corona.pos[0]+rho*cos(phi);
  pos[1] = Corona.pos[1]+rho*sin(phi);
  pos[2] = Corona.pos[2]+r*sin(theta);
}

void CoronaFindToroidalBorder(CoronaRegion Corona,double& tpre,
                              double pos[],double direct[])
{
  double major,minor,t,LocalPos[3];

  major = Corona.radius;
  minor = CoronaToroidalTubeRadius(Corona);
  t = tpre;
  VecDiff(LocalPos,pos,Corona.pos);

  if (major>0 && minor>0)
  {
    TorusCoinc(LocalPos,direct,t,minor,minor,major,0);
  }

  if (t<tpre)
  {
    tpre = t;
  }
}


bool CoronaInsideToroidal(CoronaRegion Corona,double pos[],double gap)
{
  double LocalPos[3],major,minor,eps,rho,dr,tube;

  VecDiff(LocalPos,pos,Corona.pos);
  major = Corona.radius;
  minor = CoronaToroidalTubeRadius(Corona);
  eps = CoronaBoundaryEps(major,minor,gap);
  tube = minor+eps;
  rho = sqrt(LocalPos[0]*LocalPos[0]+LocalPos[1]*LocalPos[1]);
  dr = rho-major;

  return (major>0 && minor>0 &&
          dr*dr+LocalPos[2]*LocalPos[2]<=tube*tube);
}
