// Slab corona geometry helpers.

#ifndef corona_slab_h
#define corona_slab_h

#include "model-stokes-v1.2.h"

void CoronaApplySlabDefaults(CoronaRegion& Corona,XRBSeedDisk& SeedDisk,
                             bool readCoronaSource,bool readSeedSource,
                             bool readSeedAngular,bool readDiskReflection);

void CoronaSampleSlab(CoronaRegion Corona,double r1,double r2,
                      double r3,double pos[]);

void CoronaFindSlabBorder(CoronaRegion Corona,double& tpre,
                          double pos[],double direct[]);

bool CoronaInsideSlab(CoronaRegion Corona,double pos[],double gap);

#endif
