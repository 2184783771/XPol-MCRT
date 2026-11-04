// Lamp-post corona geometry helpers.

#ifndef corona_lamppost_h
#define corona_lamppost_h

#include "model-stokes-v1.2.h"

void CoronaApplyLampPostDefaults(CoronaRegion& Corona,XRBSeedDisk& SeedDisk,
                                 bool readCoronaSource,bool readSeedSource,
                                 bool readSeedAngular,bool readDiskReflection);

void CoronaSampleLampPost(CoronaRegion Corona,double r1,double r2,
                          double r3,double pos[]);

void CoronaFindLampPostBorder(CoronaRegion Corona,double& tpre,
                              double pos[],double direct[]);

bool CoronaInsideLampPost(CoronaRegion Corona,double pos[],double gap);

#endif
