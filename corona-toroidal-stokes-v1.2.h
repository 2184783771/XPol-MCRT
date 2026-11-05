// Toroidal X-ray corona geometry.
#ifndef CORONA_TOROIDAL_H
#define CORONA_TOROIDAL_H

#include "model-stokes-v1.2.h"

void CoronaApplyToroidalDefaults(CoronaRegion& Corona,XRBSeedDisk& SeedDisk,
                                 bool readCoronaSource,bool readSeedSource,
                                 bool readSeedAngular,bool readDiskReflection);
void CoronaSampleToroidal(CoronaRegion Corona,double r1,double r2,
                          double r3,double pos[]);
void CoronaFindToroidalBorder(CoronaRegion Corona,double& tpre,
                              double pos[],double direct[]);
bool CoronaInsideToroidal(CoronaRegion Corona,double pos[],double gap);

#endif
