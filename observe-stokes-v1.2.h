// Header file listing the routines for the virtual detector system.
// The object 'observe' manages a spherical web of detectors. It
// contains routines to register photons escaping from the model
// region and to save the acquired polarimetric data of each detector
// to a disk.
//
// STOKES, version 1.2, March 2014


#ifndef observe_h
#define observe_h

#include "detector-stokes-v1.2.h"


// type definition of the detector system
class Monitor
{

public:


// initialize detectors by setting all channels to zero
Monitor();


// register photon by the detector in its flight direction
// and save its wavelength and Stokes vector
void RegisterPhoton(Photon p,Model M);


// save the data of all detectors to a disk
void SaveResults(long int Photons,Model M);


private:


// type definition of a ring of detectors for a given theta angle
struct Column
{
  Screen Row[ScreenNumPhiMax];
};


// type definition of a pointer for Column
typedef Column* ColumnPtr;


// declare a line of detectors as a dynamical array
ColumnPtr Col;
ColumnPtr DirectCol;
ColumnPtr ReflectedCol;
ColumnPtr CoronaScatteredCol;


// register photon into a selected detector set
void RegisterPhotonToColumn(Photon p,Model M,ColumnPtr Data);


// save a selected detector set to disk
void SaveDetectorSet(long int Photons,Model M,ColumnPtr Data,
                     const char suffix[]);

// save one band-integrated spatial Stokes map to disk
void SaveMapSet(long int Photons,Model M,ColumnPtr Data,
                const char suffix[]);

void SaveBandSummary(long int Photons,Model M,ColumnPtr Data,
                     const char suffix[]);


// angular difference between two neighboring viewing angles in phi
// and in theta
double PhiStep,ThetaStep;


};


#endif
