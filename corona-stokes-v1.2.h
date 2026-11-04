// Corona geometry and source sampling helpers.
//
// The model-specific geometry is split into one source file per corona
// type, while this header exposes a small dispatch layer used by Photon.

#ifndef corona_h
#define corona_h

#include "model-stokes-v1.2.h"
#include "corona-slab-stokes-v1.2.h"
#include "corona-sphere-stokes-v1.2.h"
#include "corona-lamppost-stokes-v1.2.h"

void CoronaApplySourceDefaults(CoronaRegion& Corona,XRBSeedDisk& SeedDisk,
                               bool readCoronaSource,bool readSeedSource,
                               bool readSeedAngular,bool readDiskReflection);

void CoronaSamplePosition(CoronaRegion Corona,double r1,double r2,
                          double r3,double pos[]);

double CoronaSampleEnergy(CoronaRegion Corona,double r);

void CoronaFindNextBorder(CoronaRegion Corona,double& tpre,
                          double pos[],double direct[]);

bool CoronaInside(CoronaRegion Corona,double pos[],double gap);

double CoronaBoundaryEps(double radius,double height,double gap);

#endif
