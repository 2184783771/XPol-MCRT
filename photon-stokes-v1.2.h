// Header file listing all routines concerning a model photon. The
// object 'photon follows a photon's life from its creation in some
// emission region of the model space, through scattering processes,
// until either its absorption or detection. All necessary routines
// and parameters are defined here, including the procedures for
// electron scattering
// grain. Some features are marked with a symbol (*) denoting that
// they are not used by the program at its current state but that they
// are reserved for future purposes or diagnostic reasons.
//
// STOKES, version 1.2, March 2014


#ifndef photon_h
#define photon_h


#include "model-stokes-v1.2.h"


// type definition for a photon
class Photon
{


public:


  // set the photon parameters position, flying direction, Stokes
  // vector, weight, and wavelength to the given values (*)
  void SetPhotonData(double setpos[],double setdirect[][3],double setI,
                     double setQ,double setU,double setV,double setW,
                     double setlambda);


  // reset photon to an initial state, located in a continuum or line
  // emission region, being unpolarized, and flying toward a random
  // direction
  void Reset(Model M);


  // reset photon to an initial state emitted by the corona
  void ResetFromCorona(Model M);

  // reset photon to an initial state emitted by the XRB disk
  void ResetFromXRBSeedDisk(Model M);


  // sample a photon energy in keV from the corona power-law spectrum
  double GetCoronaEnergy(CoronaRegion Corona,double r);


  // sample a wavelength for the photon, decide according to the
  // parameter 'line' if the photon is a broad or a narrow line photon
  double GetLineLambda(Model M,int line,double r);


  // if the scattering region Cloud is in the photon's flight
  // direction, then set 'tpre' to the distance to the cloud's closest
  // border, otherwise leave 'tpre' unchanged
  void FindNextBorder(scatterer Cloud,double& tpre,double pos[],
                      double direct[]);


  // if the corona is in the photon's flight direction, then set
  // 'tpre' to the distance to the closest corona border
  void FindNextCoronaBorder(CoronaRegion Corona,double& tpre,
                            double pos[],double direct[]);


  // search for the nearest scattering region along the photons flight
  // direction, if there is such cloud, then move the photon to the
  // clouds border
  bool ShiftToNextScatt(Model M);


  // sample a path length for the photon and shift it within the given
  // scattering region. If the photon leaves the border of the region
  // return false, otherwise return true
  bool Proceed(Model M);


  // return true if the photon is located inside the scattering region Cloud,
  // return false otherwise
  bool InsideCloud(scatterer Cloud,double gap=0);


  // return true if the photon is located inside the corona region
  bool InsideCorona(CoronaRegion Corona,double gap);


  // return true if the current scattering region is the corona
  bool InCorona();


  // return true if the photon would cross the z=0 XRB disk before
  // reaching infinity; set radius to the disk-plane crossing radius
  bool WouldHitXRBDisk(Model M,double& radius);


  // return true if the photon has crossed the XRB disk plane during
  // propagation; set radius to the first crossing radius
  bool HasHitXRBDisk(double& radius);


  // return the number of XRB disk reflections undergone by the photon
  int GetXRBDiskReflectionNum();


  // return the number of corona scatterings undergone by the photon
  int GetCoronaScatterNum();


  // return true if the photon belongs to the reflected disk component
  bool IsReflectedComponent();


  // return true if the photon belongs to the direct component
  bool IsDirectComponent();


  // return true if the photon was absorbed by the XRB disk reflection layer
  bool AbsorbedByXRBDisk();


  // reflect an escaped photon from the XRB disk if its current
  // trajectory intersects the disk before reaching the detector
  bool ReflectEscapedXRBDisk(Model M);


  // turn the photon's flight direction by the angle 'alpha' and around the
  // coordinate axis given by int
  void TurnDirection(int axis,double alpha);


  // turn the photon's polar plane by the angle 'alpha' around the axis given
  // by its flight direction
  void TurnPolarPlane(double alpha);


  // perform an electron scattering event and change all affected
  // photon parameters
  void ElectronScattering(Model M);


  // sample an electron scattering angle theta according to the
  // random number r
  double GetElecScattTheta(double r);


  // sample an electron scattering angle phi according to the
  // random number r, incident polarization PL, and the scattering
  // angle theta
  double GetElecScattPhi(double r,double PL,double theta);


  // turn the flight direction of the photon due to electron scattering by
  // the angle theta, adjust the Stokes parameters
  void ElectTurnTheta(double theta);


  // transform the photon wavelength to (ToFrame=true) or from (ToFrame=false)
  // a reference frame moving in the radial direction with the velocity vr
  void LambdaShift(bool ToFrame,int vtype,double vr,double vphi,double vtheta); 


  // adjust the Stokes parameters of the photon to the polarization plane of
  // a detector in the direction given by theta and phi
  void AdjustFrame(double phi,double theta);


  // return the photon wavelength
  double GetLambda();


  // return the X resolution for the imagery plane
  double GetResX();


  // return the Y resolution for the imagery plane
  double GetResY();


  // set photon position to 'setpos'
  void SetPosition(double setpos[]);


  // set 'getpos' to the current photon position
  void GetPosition(double getpos[]);


  // return the x-component of the flight direction
  double GetDirect1();


  // return the y-component of the flight direction
  double GetDirect2();


  // return the z-component of the flight direction
  double GetDirect3();


  // return the x-component of the photon position
  double GetPos1();


  // return the y-component of the photon position
  double GetPos2();


  // return the z-component of the photon position
  double GetPos3();


  // set the scattering region number to 'NewFlag'
  void SetScattFlag(int NewFlag);


  // return the current scattering region number
  int GetFlag();


  // return the number of scatterings the photon has undergone (*)
  int GetScattNum();


  // return the distance of the photon from the origin of the model space
  double GetDistance();


  // return the I value of the Stokes vector
  double GetI();


  // return the Q value of the Stokes vector
  double GetQ();


  // return the U value of the Stokes vector
  double GetU();


  // return the V value of the Stokes vector
  double GetV();


  // return the weight factor of the photon (*)
  double GetW();


  // return the value of the photon's absorption flag
  bool Absorbed();


  // set the value of the photon's absorption flag to true
  void SetAbsorbed();


private:


  // reflect from the z=0 XRB disk if the segment crosses the disk
  bool ReflectFromXRBDisk(Model M,double OldPos[],double NewPos[]);

  // set a new propagation frame from a flight direction
  void SetDirectionFrame(double NewDirect[]);


  // number of the current scattering region
  int ScattFlag;

  // number of scatterings the photon has undergone (*)
  int ScattNum;

  // position and flight direction
  double pos[3],direct[3][3];

  // Stokes vector, wavelength, and weight
  // (the weight factor is reserved for future purposes)
  double I,Q,U,V,lambda,W;

  // cross-section of the scattering media in the
  // current cloud
  double CS;

  // absorption flag
  bool abs;

  // diagnostic flag for crossing the XRB disk before detection
  bool XRBDiskHit;
  double XRBDiskHitRadius;
  int XRBDiskReflectionNum;
  int XRBDiskHitNum;
  int CoronaScatterNum;
  int OriginFlag;
  bool XRBDiskAbsorbed;

  // path length after changing scattering region
  double PathLength;

};


#endif
