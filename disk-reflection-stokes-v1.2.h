// XRB disk reflection helpers.
//
// This module handles intersection with the flat z=0 disk and computes
// a simple single-scattering disk reflection event. Photon owns the
// photon state and history; this module only returns the event result.

#ifndef disk_reflection_h
#define disk_reflection_h

#include "model-stokes-v1.2.h"

struct DiskReflectionEvent
{
  bool hit;
  bool reflected;
  bool absorbed;
  double radius;
  double pos[3];
  double direct[3];
  double polarAxis[3];
  double Q,U,V;
};

void InitDiskReflectionEvent(DiskReflectionEvent& Event);

bool DiskReflectionIntersect(Model M,double OldPos[],double NewPos[],
                             double& radius);

void DiskReflectionApply(Model M,double OldPos[],double NewPos[],
                         double Incoming[],double Frame1[],
                         double Frame2[],double I,double Q,double U,
                         double V,DiskReflectionEvent& Event);

#endif
