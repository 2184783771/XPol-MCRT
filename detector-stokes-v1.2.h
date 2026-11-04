// Header file listing the routines for virtual detector. The object
// 'detector' adds up the Stokes vectors of the detected photons
// according to their wavelength. It then can save the stored data to
// a disk.
//
// STOKES, version 1.2, March 2014


#ifndef detector_h
#define detector_h


#include <fstream>
#include "model-stokes-v1.2.h"
#include "photon-stokes-v1.2.h"

using namespace std;



// type definition for a detector
class Screen
{


public:


  // initialize detector by setting all channels to zero
  Screen();


  // register a photon according to its wavelength and add its
  // Stokes parameters
  // save de XY directions for the imaging plane
  void RegisterPhoton(Photon p,Model M, double x, double y);


  // save the data of the detector to the disk
  void SaveCurve(ofstream out_stream, int n, Model M);


  // return a stored Stokes value
  double GetStokes(int n,int channel,int x,int y);


private:


  // allocate detector storage for the actual model resolution
  void EnsureAllocated(Model M);


  // flatten detector indices
  int Index(int channel,int x,int y);


  // dynamically allocated Stokes channel arrays
  double *I,*Q,*U,*V;
  int AllocRes,AllocX,AllocY;


};


#endif
