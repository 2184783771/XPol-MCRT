// Header file containing the definition of the 'model' routines. The
// object 'model' stores all parameters for the model run defined in
// 'config.ini'. Basic parameters and structures of the model space,
// such as emission and scattering regions, are defined here. The
// routines of the object 'model' read in the input data from the file
// 'config.ini' and compute basic parameters of the specified run.
//
// STOKES, version 1.2, March 2014


#ifndef model_h
#define model_h

#include <iosfwd>
#include "specmath-stokes-v1.2.h"


// maximum number of spectroscopic channels of the detectors
const int ResMax = 240;


// maximum resolution of the X-axis imagery plane 
const int ResXmax = 128;


// maximum resolution of the Y-axis imagery plane 
const int ResYmax = 128;


// maximum number of detectors in the phi-direction
const int ScreenNumPhiMax = 40;


// maximum number of detectors in the theta-direction
const int ScreenNumThetaMax = 80;


// maximum number of scattering regions
const int ScattMax = 30;


// special scattering flag for the optional corona region
const int CoronaFlag = ScattMax + 1;


// precision for the computation of the phase functions (in degrees)
const double PFPrecision = 1.0;


// precision for the sampling of the line photons (in Angstroem)
const double LinePrecision = 0.01;


// Thomson cross section in ccm/parsec
const double TCS = 3.08633e18*6.6524e-25;


// type of a scatterer (torus, cylinder, conical, sphere) - shape,
// pos, a, b, and theta describe the geometry of the region. Density is
// the electron number density and vr denotes a radial velocity.
struct scatterer
{
  int shape,vtype;
  double pos[3];
  double R,a,b,c,theta;
  double vr,vphi,vtheta;
  double Density;
};


// type of a continuum emission region (disk-shaped), spectral index
// alpha for the intensity spectrum I(nu) ~ nu^(-alpha)
struct CTRegion
{
  int vtype;
  double R,a,b,c,alpha,vr,vphi,vtheta;
};


// type of a broad line emission region (toroidal)
struct BLRegion
{
  int vtype;
  double R,a,b,lambda,width,rli,vr,vphi,vtheta,norm,low;
};


// type of a narrow line emission region (conical)
struct NLRegion
{
  int vtype;
  double a,b,theta,lambda,width,rli,vr,vphi,vtheta,norm,low;
};


// type of a corona region reserved for X-ray corona extensions.
// system: 0 X-ray binary.
// type: 0 slab, 1 toroidal, 2 lamp-post.
struct CoronaRegion
{
  int system,type;
  bool source,compton;
  double pos[3];
  double radius,height;
  double tau,density,kTe;
  double gamma,emin,emax;
};

// type of an X-ray binary disk seed photon source.
// Photons are emitted from a flat disk in the z=0 plane.
struct XRBSeedDisk
{
  int spectrum,angular,polarization,reflectAngular;
  bool source,reflection;
  double rin,rout,reflectRin,reflectRout,albedo;
  double gamma,tin,emin,emax;
};


// type definition of the object model
class Model
{


public:


  // initialize the Model parameters by setting them to the
  // default values
  Model();


  // read in model parameters from 'config.ini', if the initialization
  // happens correctly return 'true'.
  bool Initialize();


  // return the number of detectors in the phi-direction
  int ScreenNumPhi();


  // return the number of detectors in the theta-direction
  int ScreenNumTheta();


  // return the number of spectroscopic channels
  int res();


  // return the number of photons to sample
  long int PhotonNum();


  // return the number of scattering regions defined
  int GetScattNum();


  // return parameters of the scattering region with index i
  scatterer GetCloud(int i);


  // return the maximum expansion of all scattering regions
  double GetMaxExp();


  // return the X axis resolution
  double GetResX();


  // return the Y axis resolution
  double GetResY();


  // return the maximum radial velocity of all scattering regions
  double GetVMax();


  // return the minimum photon wavelength defined
  double LambdaMin();


  // return the maximum photon wavelength defined
  double LambdaMax();


  // return the photon wavelength range covered
  double LambdaRange();

  // return the random seed configured for this run
  int GetRandomSeed();

  // return true if the random seed was explicitly read from config.ini
  bool HasFixedRandomSeed();

  // return the simplified spin case used for the ISCO-truncation test
  int GetSpinCase();

  // return true if SpinCase was explicitly read from config.ini
  bool HasSpinCase();

  // return the ISCO radius used by SpinCase in units of Rg
  double GetSpinISCO();

  // return the projected system-axis position angle on the sky, measured
  // in degrees from celestial north toward east
  double GetSkyAxisPA();

  bool HasInclination();
  double GetInclinationDeg();
  double GetSourceRA();
  double GetSourceDec();
  double GetSourceDistanceKpc();
  double GetBlackHoleMassSolar();


  // set the output file name to 'name'
  void GetOutputFile(char name[300]);


  // return the minimum gab between scattering clouds
  double GetGap();


  // return 'true' if the model space is defined symmetrically with
  // respect to the xy-plane and 'false' otherwise
  bool GetPlaneSymmetry();


  // return 'true' if the model space is defined symmetrically with
  // respect to the z-axis and 'false' otherwise
  bool GetAxisSymmetry();


  // return the parameters of the continuum emission region
  CTRegion GetCTR();


  // return the parameters of the broad line emission region
  BLRegion GetBLR();


  // return the parameters of the narrow line emission region
  NLRegion GetNLR();


  // return 'true' if a corona region is defined in config.ini
  bool HasCorona();


  // return 'true' if photons should be emitted from the corona
  bool HasCoronaSource();

  // return 'true' if XRB seed photons should be emitted from the disk
  bool HasXRBSeedSource();

  // return 'true' if circular polarization V should be saved
  bool ShouldSaveV();

  // return 'true' if direct/reflected component files should be saved
  bool ShouldSaveComponents();

  // return 'true' if full wavelength-resolved .dat files should be saved
  bool ShouldSaveStokesFiles();

  // return 'true' if band-integrated spatial map CSV files should be saved
  bool ShouldSaveMap();

  // return the detector theta index used for spatial map CSV output
  int GetMapThetaIndex();

  // return the energy band used for spatial map CSV output
  double GetMapEnergyMin();
  double GetMapEnergyMax();


  // return the parameters of the corona region
  CoronaRegion GetCorona();

  // return the parameters of the XRB disk seed source
  XRBSeedDisk GetXRBSeedDisk();


private:

  bool ParseGeneralOption(std::istream& input,char token[200]);
  bool ParseLegacyOption(std::istream& input,char token[200]);
  bool ParseCoronaOption(std::istream& input,char token[200],
                         bool& readCoronaSource);
  bool ParseDiskOption(std::istream& input,char token[200],
                       bool& readSeedSource,bool& readSeedAngular,
                       bool& readDiskReflection,bool& readDiskReflectRin,
                       bool& readDiskReflectRout);
  bool ParseOutputOption(std::istream& input,char token[200]);


  // parameters of the continuum emission region
  CTRegion CTR;


  // parameters of the broad emission line region
  BLRegion BLR;


  // parameters of the narrow emission line region
  NLRegion NLR;


  // parameters of the corona region
  CoronaRegion Corona;


  // set to true when corona parameters are read from config.ini
  bool UseCorona;

  // parameters of the XRB disk seed photon source
  XRBSeedDisk SeedDisk;


  // parameters of the scattering regions
  scatterer Cloud[ScattMax];


  // number of scattering regions
  int EndScattList;


  // number of detectors in the phi-direction
  int SNumPhi;


  // number of detectors in the theta direction
  int SNumTheta;


  // number of spectroscopic channels for the detectors
  int resolut;
 

  // minimum and maximum wavelength (in Angstroem) of photons sampled
  double LamMin,LamMax;
  

  // resolution of X and Y axis
  double ResX,ResY;


  // total number of photons sampled
  long int PNum;

  // Random-number control.  If RandomSeedFixed is false, the program keeps the
  // historical time-based seed initialized by the photon module.
  int RandomSeed;
  bool RandomSeedFixed;

  // Simplified spin control: SpinCase 0 -> a*=0, R_ISCO=6 Rg;
  // SpinCase 1 -> a*=0.998, R_ISCO=1.24 Rg.
  int SpinCase;
  bool SpinCaseEnabled;
  double SpinISCO;

  // Projected disk/spin/jet axis position angle on the sky, in degrees
  // east of celestial north.  This rotates only the detector Q/U basis.
  double SkyAxisPA;
  double InclinationDeg;
  bool InclinationEnabled;
  double SourceRA,SourceDec,SourceDistanceKpc,BlackHoleMassSolar;


  // maximum expansion among all scattering regions
  double MaxExp;


  // maximum radial velocity among all scattering regions
  double VMax;


  // Minimum gap between scattering regions.This parameter
  // is required for computational reasons. The object Photon
  // distinguishes between photon positions inside and outside
  // scattering clouds. A photon that is leaving a scattering cloud is
  // re-located outside of it, at the distance GetGap to the (smoothly
  // defined) cloud surface. This is done to ensure that the photon is
  // safely recognized as being outside the cloud.
  double Gap;


  // Stem name of the output files
  char outputfile[80];

  // Parameter describing if the model space defined is symmetric with
  // respect to the xy-plane. If it is set to 'true' the results at the
  // inclinations i and 180 deg - i are combined.
  bool planesym;

  // Parameter describing if the model space defined is symmetric with
  // respect to the z-axis. If it is set to 'true' the results at all
  // azimuthal angles are combined.
  bool axissym;

  // Output controls.  Full STOKES compatibility keeps both enabled by default,
  // while XRB Stokes-map diagnostics can disable V to reduce disk I/O.
  bool SaveV;
  bool SaveComponents;
  bool SaveStokesFiles;
  bool SaveMap;
  int MapThetaIndex;
  double MapEnergyMin,MapEnergyMax;

};


#endif
