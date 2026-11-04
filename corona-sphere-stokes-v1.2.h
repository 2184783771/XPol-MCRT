// Compact spherical helpers for lamp-post sources, plus the toroidal
// XRB corona geometry used by CoronaType toroidal.

#ifndef corona_sphere_h
#define corona_sphere_h

#include "model-stokes-v1.2.h"

void CoronaApplySphereDefaults(CoronaRegion& Corona,XRBSeedDisk& SeedDisk,
                               bool readCoronaSource,bool readSeedSource,
                               bool readSeedAngular,bool readDiskReflection);

void CoronaApplyToroidalDefaults(CoronaRegion& Corona,XRBSeedDisk& SeedDisk,
                                 bool readCoronaSource,bool readSeedSource,
                                 bool readSeedAngular,bool readDiskReflection);

void CoronaSampleSphere(CoronaRegion Corona,double r1,double r2,
                        double r3,double pos[]);

void CoronaSampleToroidal(CoronaRegion Corona,double r1,double r2,
                          double r3,double pos[]);

void CoronaFindSphereBorder(CoronaRegion Corona,double& tpre,
                            double pos[],double direct[]);

void CoronaFindToroidalBorder(CoronaRegion Corona,double& tpre,
                              double pos[],double direct[]);

bool CoronaInsideSphere(CoronaRegion Corona,double pos[],double gap);

bool CoronaInsideToroidal(CoronaRegion Corona,double pos[],double gap);

#endif
