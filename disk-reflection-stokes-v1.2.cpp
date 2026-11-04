#include <math.h>
#include "disk-reflection-stokes-v1.2.h"
#include "disk-reflection-polarization-stokes-v1.2.h"
#include "randomc-stokes-v1.2.h"
#include "specmath-stokes-v1.2.h"

extern CRandomMersenne RanGen;


void InitDiskReflectionEvent(DiskReflectionEvent& Event)
{
  int i;

  Event.hit = false;
  Event.reflected = false;
  Event.absorbed = false;
  Event.radius = 0;
  Event.Q = 0;
  Event.U = 0;
  Event.V = 0;
  for(i=0; i<=2; i++)
  {
    Event.pos[i] = 0;
    Event.direct[i] = 0;
    Event.polarAxis[i] = 0;
  }
}


bool DiskReflectionIntersect(Model M,double OldPos[],double NewPos[],
                             double& radius)
{
  double f,x,y;
  XRBSeedDisk Disk;

  radius = 0;
  Disk = M.GetXRBSeedDisk();

  if (Disk.reflectRout<=Disk.reflectRin || Disk.reflectRout<=0)
  {
    return false;
  }

  // Count upper-side photons crossing downward through the disk plane.
  if (OldPos[2]<=0 || NewPos[2]>0)
  {
    return false;
  }

  f = OldPos[2]/(OldPos[2]-NewPos[2]);
  if (f<0 || f>1)
  {
    return false;
  }

  x = OldPos[0]+f*(NewPos[0]-OldPos[0]);
  y = OldPos[1]+f*(NewPos[1]-OldPos[1]);
  radius = sqrt(x*x+y*y);

  return (radius>=Disk.reflectRin && radius<=Disk.reflectRout);
}


void DiskReflectionApply(Model M,double OldPos[],double NewPos[],
                         double Incoming[],double Frame1[],
                         double Frame2[],double I,double Q,double U,
                         double V,DiskReflectionEvent& Event)
{
  double f,x,y,mu,phi,st;
  DiskReflectionPolarization Polarization;
  XRBSeedDisk Disk;

  InitDiskReflectionEvent(Event);
  Disk = M.GetXRBSeedDisk();
  if (!DiskReflectionIntersect(M,OldPos,NewPos,Event.radius))
  {
    return;
  }

  Event.hit = true;
  f = OldPos[2]/(OldPos[2]-NewPos[2]);
  x = OldPos[0]+f*(NewPos[0]-OldPos[0]);
  y = OldPos[1]+f*(NewPos[1]-OldPos[1]);
  Event.pos[0] = x;
  Event.pos[1] = y;
  Event.pos[2] = M.GetGap();

  if (!Disk.reflection)
  {
    return;
  }
  if (RanGen.Random()>Disk.albedo)
  {
    Event.absorbed = true;
    return;
  }

  if (Disk.reflectAngular==1)
  {
    mu = sqrt(RanGen.Random());
    phi = 2*pi*RanGen.Random();
    st = sqrt(1-mu*mu);
    Event.direct[0] = st*cos(phi);
    Event.direct[1] = st*sin(phi);
    Event.direct[2] = mu;
  }
  else
  {
    Event.direct[0] = Incoming[0];
    Event.direct[1] = Incoming[1];
    Event.direct[2] = -Incoming[2];
  }

  if (!DiskReflectionPolarize(Incoming,Event.direct,Frame1,Frame2,
                              I,Q,U,V,Polarization))
  {
    Event.absorbed = true;
    return;
  }
  Event.polarAxis[0] = Polarization.axis[0];
  Event.polarAxis[1] = Polarization.axis[1];
  Event.polarAxis[2] = Polarization.axis[2];
  Event.Q = Polarization.Q;
  Event.U = Polarization.U;
  Event.V = Polarization.V;
  Event.reflected = true;
}
