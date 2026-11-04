// Implementation file for the 'model' routines. The object 'model'
// reads in and stores all parameters for the model run defined in
// 'config.ini'. Basic parameters and structures of the model space,
// such as emission and scattering regions, are defined here. The
// routines of the object 'model' read in the input data from the file
// 'config.ini' and compute basic parameters of the specified model run.
//
// STOKES, version 1.2, March 2014


#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "model-stokes-v1.2.h"
#include "corona-stokes-v1.2.h"

using namespace std;



// Initialize the Model parameters by defining an empty model
// space. Assume a flat intensity spectrum over the whole possible
// wavelength range of 1000 Ang to 10000 Ang coming from a point-like
// continuum emission region. The (not used) line emission regions are
// set to a cylindrical broad line region and to  a double conical
// narrow line region for the H-alpha line at 6563 Ang. The (not used)
// theta direction ten viewing angles are defined, a spectroscopic
// view with 45 channels is declared. A million photons are sampled.
Model::Model()
{
  int i,j;

  PNum = 1000000;
  RandomSeed = 0;
  RandomSeedFixed = false;
  SpinCase = -1;
  SpinCaseEnabled = false;
  SpinISCO = 0;
  SkyAxisPA = 0;
  InclinationDeg = 60.0;
  InclinationEnabled = false;
  SourceRA = 0.0;
  SourceDec = 0.0;
  SourceDistanceKpc = 0.0;
  BlackHoleMassSolar = 0.0;
  LamMin = 1000;
  LamMax = 10000;
  ResX = 1;
  ResY = 1;
  SaveV = true;
  SaveComponents = true;
  SaveStokesFiles = true;
  SaveMap = false;
  MapThetaIndex = 0;
  MapEnergyMin = 2.0;
  MapEnergyMax = 8.0;
  CTR.R = 0;
  CTR.a = 0;
  CTR.b = 0;
  CTR.c = 0;
  CTR.alpha = 2;
  CTR.vtype = 0;
  CTR.vr = 0;
  CTR.vphi = 0;
  CTR.vtheta = 0;
  BLR.R = 0.0070;
  BLR.a = 0.0050;
  BLR.b = 0.0050;
  BLR.lambda = 6563;
  BLR.width = 150;
  BLR.rli = 0;
  BLR.vtype = 0;
  BLR.vr = 0;
  BLR.vphi = 0;
  BLR.vtheta = 0;
  NLR.a = 10;
  NLR.b = 100;
  NLR.theta = pi/180*25;
  NLR.lambda = 6563;
  NLR.width = 25;
  NLR.rli = 0;
  NLR.vtype = 0;
  NLR.vr = 0;
  NLR.vphi = 0;
  NLR.vtheta = 0;
  UseCorona = false;
  Corona.system = 0;
  Corona.type = 1;
  Corona.source = false;
  Corona.compton = false;
  for(j=0; j<=2; j++) Corona.pos[j] = 0;
  Corona.radius = 0;
  Corona.height = 0;
  Corona.tau = 0;
  Corona.density = 0;
  Corona.kTe = 0;
  Corona.gamma = 2;
  Corona.emin = 1;
  Corona.emax = 100;
  SeedDisk.spectrum = 0;
  SeedDisk.angular = 0;
  SeedDisk.polarization = 0;
  SeedDisk.reflectAngular = 0;
  SeedDisk.source = false;
  SeedDisk.reflection = false;
  SeedDisk.rin = 0;
  SeedDisk.rout = 10;
  SeedDisk.reflectRin = 0;
  SeedDisk.reflectRout = 0;
  SeedDisk.albedo = 1;
  SeedDisk.gamma = 2;
  SeedDisk.tin = 1;
  SeedDisk.emin = 1;
  SeedDisk.emax = 10;
  for(i=0; i<ScattMax; i++)
  {
    Cloud[i].shape = 0;
    for(j=0; j<=2; j++) Cloud[i].pos[j] = 0;
    Cloud[i].R = 0;
    Cloud[i].a = 0;
    Cloud[i].b = 0;
    Cloud[i].c = 0;
    Cloud[i].theta = 0;
    Cloud[i].vtype = 0;
    Cloud[i].vr = 0;
    Cloud[i].vphi = 0;
    Cloud[i].vtheta = 0;
    Cloud[i].Density = 0;
  }
  EndScattList = 0;
  Gap = 1.0e-9;
  SNumPhi = 1;
  SNumTheta = 10;
  planesym = false;
  axissym = false;
  resolut = 45;
  strcpy(outputfile,"empty-space1");

  // set output format for numbers
  cout.setf(ios::fixed);
  cout.setf(ios::showpoint);
  cout.precision(5);

  // welcome screen
  cout << "\n\n";
  cout << "========================================\n\n";
  cout << "               XPol-MCRT\n\n";
  cout << " X-ray Polarization Monte Carlo Radiative Transfer\n\n";
  cout << " Based on STOKES v1.2 by Rene Goosmann\n";
  cout << " and Frederic Marin\n";
  cout << "========================================\n\n";
}



// read in model parameters from 'config.ini', if the initialization
// happens correctly return 'true'.
bool Model::Initialize()
{
  bool accept;
  bool ReadCoronaSource,ReadSeedSource,ReadSeedAngular;
  bool ReadDiskReflection,ReadDiskReflectRin,ReadDiskReflectRout;
  int i;
  double R,a,b,d;
  char str[200];
  fstream in_stream;

  in_stream.open("config.ini",ios::in);
#ifdef XPOL_SOURCE_DIR
  if (in_stream.fail())
  {
    in_stream.clear();
    const std::string sourceConfig = std::string(XPOL_SOURCE_DIR) + "/config.ini";
    in_stream.open(sourceConfig.c_str(),ios::in);
  }
#endif
  if (in_stream.fail()!=0)
  {
    cout << "Error - unable to open file 'config.ini'!\n";
    cout << "Program aborted...\n\n\n";
    return false;
  }
  else
  {

    ReadCoronaSource = false;
    ReadSeedSource = false;
    ReadSeedAngular = false;
    ReadDiskReflection = false;
    ReadDiskReflectRin = false;
    ReadDiskReflectRout = false;
    
    do
    {
      in_stream >> str;

      // in case of reaching end of file
      if (in_stream.eof())
      {
        strcpy(str,"end");
      }

      // first assume 'str' is not a key word
      accept = false;

      // in case 'str' is a comment go to next line
      if (str[0]=='#')
      {
        in_stream.getline(str,200);
	accept = true;
      }

      // else read in the key word and associated data
      else
      {
        accept = ParseGeneralOption(in_stream,str);
        if (ParseLegacyOption(in_stream,str)) accept = true;
        if (ParseCoronaOption(in_stream,str,ReadCoronaSource)) accept = true;
        if (ParseDiskOption(in_stream,str,ReadSeedSource,ReadSeedAngular,
                            ReadDiskReflection,ReadDiskReflectRin,
                            ReadDiskReflectRout)) accept = true;
        if (ParseOutputOption(in_stream,str)) accept = true;
      }
    }
    while((accept)&&(strcmp(str,"end")!=0));
    in_stream.close();

    // If data read in correctly compute some global parameters
    if (strcmp(str,"end")==0)
    {
      // Apply XRB corona model defaults unless the user explicitly
      // chose a photon source or seed angular distribution.
      if (UseCorona && Corona.system==0)
      {
        CoronaApplySourceDefaults(Corona,SeedDisk,ReadCoronaSource,
                                  ReadSeedSource,ReadSeedAngular,
                                  ReadDiskReflection);
      }
      if (!ReadDiskReflectRin)
      {
        SeedDisk.reflectRin = SeedDisk.rin;
      }
      if (!ReadDiskReflectRout)
      {
        SeedDisk.reflectRout = SeedDisk.rout;
      }
      if (SpinCaseEnabled)
      {
        SeedDisk.rin = SpinISCO;
        SeedDisk.reflectRin = SpinISCO;
      }
      // In the case of a mapping routine, the model cannot be planar symmetric
      // and the value of PlaneSym is set to "no". The user is then notified.
      if((ResX > 1 || ResY > 1) && planesym == 1)
      {
	planesym = 0;
	cout << " WARNING : for imaging purpose, the model cannot be planar symetric" << endl;
	cout << "           PlaneSym is automatically set to 'no'" << endl;
      }

      // Determine the maximum expansion among all scatterers and their
      // maximum radial velocity. This provides an upper limit of the
      // overall size and dynamics of the model space as a whole.
      MaxExp = 0;
      VMax = 0;
      for(i=0; i<EndScattList; i++)
      {
	switch(Cloud[i].shape)
	{
	  case 0:
	    a = Cloud[i].a;
	    b = Cloud[i].b;
	    R = Cloud[i].R;
	    d = 2*(sqrt((R+b)*(R+b)+a*a))+Gap;
	    if (d > MaxExp) MaxExp = d;
	    break;
          case 1:
	    a = Cloud[i].a;
	    b = Cloud[i].b;
  	    R = Cloud[i].R;
  	    d = 2*(sqrt((R+b)*(R+b)+a*a))+Gap;
  	    if (d > MaxExp) MaxExp = d;
  	    break;
          case 2:  d = 2*Cloud[i].b+Gap;
    	    if (d > MaxExp) MaxExp = d;
  	    break;
 	  case 3:  d = 2*(sqrt(SProd(Cloud[i].pos,Cloud[i].pos))+Cloud[i].a)+Gap;
            if (d > MaxExp) MaxExp = d;
   	    break;
          default: d = 2*Cloud[i].b+Gap;
    	    if (d > MaxExp) MaxExp = d;
        }
        if (Cloud[i].vr>VMax)
        {
  	  VMax = sqrt(Cloud[i].vr*Cloud[i].vr+
                      Cloud[i].vphi*Cloud[i].vphi+
                      Cloud[i].vtheta*Cloud[i].vtheta);
        }
      }
      if (UseCorona)
      {
        switch(Corona.type)
        {
          case 0:
            d = 2*(sqrt(SProd(Corona.pos,Corona.pos))+
                   sqrt(Corona.radius*Corona.radius+
                        0.25*Corona.height*Corona.height))+Gap;
            break;
          default:
            d = 2*(sqrt(SProd(Corona.pos,Corona.pos))+Corona.radius)+Gap;
            break;
        }
        if (d > MaxExp) MaxExp = d;
      }
      if (MaxExp==0) 		// In the case there is an offset on the emission region for a normalisation spectra
      {	
	  MaxExp = 1e6;
      }

      // Determine parts of the sampling integrals for the line photons
      // in order to shorten the computations in the main loop
      BLR.low = 0.5*log(pow(LamMin-BLR.lambda,2)+pow(BLR.width,2))+
                BLR.lambda/BLR.width*atang(LamMin-BLR.lambda,BLR.width);
      BLR.norm = 0.5*log(pow(LamMax-BLR.lambda,2)+pow(BLR.width,2))+
                 BLR.lambda/BLR.width*atang(LamMax-BLR.lambda,BLR.width)-
                 BLR.low;
      NLR.low = 0.5*log(pow(LamMin-NLR.lambda,2)+pow(NLR.width,2))+
                NLR.lambda/NLR.width*atang(LamMin-NLR.lambda,NLR.width);
      NLR.norm = 0.5*log(pow(LamMax-NLR.lambda,2)+pow(NLR.width,2))+
                 NLR.lambda/NLR.width*atang(LamMax-NLR.lambda,NLR.width)-
                 NLR.low;

      if (UseCorona)
      {
        if (Corona.radius<=0)
        {
          cout << "Error - CoronaRadius must be greater than 0.\n";
          cout << "Program aborted...\n\n\n";
          return false;
        }
        if (Corona.type==0 && Corona.height<=0)
        {
          cout << "Error - CoronaHeight must be greater than 0 for CoronaType slab.\n";
          cout << "Add for example: CoronaHeight 2\n";
          cout << "Program aborted...\n\n\n";
          return false;
        }
        if (Corona.tau<=0 && Corona.density<=0)
        {
          cout << "Error - corona scattering needs CoronaTau or CoronaDensity greater than 0.\n";
          cout << "Add for example: CoronaTau 1\n";
          cout << "Program aborted...\n\n\n";
          return false;
        }
        if ((Corona.source || SeedDisk.source) &&
            (Corona.emin<=0 || Corona.emax<=0))
        {
          cout << "Error - CoronaEnergyRange values must be greater than 0.\n";
          cout << "Program aborted...\n\n\n";
          return false;
        }
        if (Corona.compton && Corona.kTe<=0)
        {
          cout << "Error - CoronaTemp must be greater than 0 when CoronaCompton is yes.\n";
          cout << "Add for example: CoronaTemp 50\n";
          cout << "Program aborted...\n\n\n";
          return false;
        }
      }
      if (SeedDisk.source)
      {
        if (SeedDisk.rout<=SeedDisk.rin || SeedDisk.rout<=0)
        {
          cout << "Error - XRBSeedRout must be greater than XRBSeedRin and greater than 0.\n";
          cout << "Program aborted...\n\n\n";
          return false;
        }
        if (SeedDisk.rin<0)
        {
          cout << "Error - XRBSeedRin cannot be negative.\n";
          cout << "Program aborted...\n\n\n";
          return false;
        }
        if (SeedDisk.emin<=0 || SeedDisk.emax<=0)
        {
          cout << "Error - XRBSeedEnergyRange values must be greater than 0.\n";
          cout << "Program aborted...\n\n\n";
          return false;
        }
        if (SeedDisk.spectrum==1 && SeedDisk.tin<=0)
        {
          cout << "Error - XRBSeedTin must be greater than 0 for XRBSeedSpectrum blackbody.\n";
          cout << "Program aborted...\n\n\n";
          return false;
        }
      }
      if (SpinCaseEnabled)
      {
        if (SeedDisk.source && SeedDisk.rout<=SpinISCO)
        {
          cout << "Error - XRBSeedRout must be greater than the SpinCase ISCO radius.\n";
          cout << "Increase XRBSeedRout above " << SpinISCO << ".\n";
          cout << "Program aborted...\n\n\n";
          return false;
        }
        if (SeedDisk.reflection && SeedDisk.reflectRout<=SpinISCO)
        {
          cout << "Error - DiskReflectRout must be greater than the SpinCase ISCO radius.\n";
          cout << "Increase DiskReflectRout above " << SpinISCO << ".\n";
          cout << "Program aborted...\n\n\n";
          return false;
        }
      }
      if (SeedDisk.reflection)
      {
        if (SeedDisk.reflectRout<=SeedDisk.reflectRin ||
            SeedDisk.reflectRout<=0)
        {
          cout << "Error - DiskReflectRout must be greater than DiskReflectRin and greater than 0.\n";
          cout << "Program aborted...\n\n\n";
          return false;
        }
        if (SeedDisk.reflectRin<0)
        {
          cout << "Error - DiskReflectRin cannot be negative.\n";
          cout << "Program aborted...\n\n\n";
          return false;
        }
        if (SeedDisk.albedo<0 || SeedDisk.albedo>1)
        {
          cout << "Error - DiskReflectAlbedo must be between 0 and 1.\n";
          cout << "Program aborted...\n\n\n";
          return false;
        }
      }
      if (resolut<1 || resolut>ResMax)
      {
        cout << "Error - SpectRes must be between 1 and "
             << ResMax << ".\n";
        cout << "Program aborted...\n\n\n";
        return false;
      }
      if (SNumTheta<1 || SNumTheta>ScreenNumThetaMax)
      {
        cout << "Error - ThetaViewAng must be between 1 and "
             << ScreenNumThetaMax << ".\n";
        cout << "Program aborted...\n\n\n";
        return false;
      }
      if (SNumPhi<1 || SNumPhi>ScreenNumPhiMax)
      {
        cout << "Error - PhiViewAng must be between 1 and "
             << ScreenNumPhiMax << ".\n";
        cout << "Program aborted...\n\n\n";
        return false;
      }
      if (InclinationEnabled && (InclinationDeg<0.0 || InclinationDeg>180.0))
      {
        cout << "Error - InclinationDeg must be between 0 and 180 degrees.\n";
        cout << "Program aborted...\n\n\n";
        return false;
      }
      if (InclinationEnabled && planesym && InclinationDeg>90.0)
      {
        cout << "Error - InclinationDeg must be at most 90 when PlaneSym is yes.\n";
        cout << "Program aborted...\n\n\n";
        return false;
      }
      if (SourceRA<0.0 || SourceRA>=360.0 || SourceDec<-90.0 || SourceDec>90.0 ||
          SourceDistanceKpc<0.0 || BlackHoleMassSolar<0.0)
      {
        cout << "Error - invalid source coordinates, distance, or black-hole mass.\n";
        cout << "Program aborted...\n\n\n";
        return false;
      }
      if (ResX<1 || ResX>ResXmax || ResY<1 || ResY>ResYmax)
      {
        cout << "Error - ResX and ResY must be between 1 and their compiled limits.\n";
        cout << "Program aborted...\n\n\n";
        return false;
      }

      cout << "Input parameters properly read in\n";
      cout << "Output file name set to '" << outputfile << "'\n";
      if (UseCorona)
      {
        cout << "Corona parameters enabled\n";
      }
      if (SeedDisk.source)
      {
        cout << "XRB disk seed photon source enabled\n";
      }
      if (SeedDisk.reflection)
      {
        cout << "XRB disk reflection enabled\n";
      }
      if (SpinCaseEnabled)
      {
        cout << "SpinCase " << SpinCase << " enabled: disk inner radius set to R_ISCO="
             << SpinISCO << " Rg\n";
      }
      if (UseCorona && Corona.system==0 && Corona.type==1 && Corona.radius>0)
      {
        cout << "Toroidal+disk geometry enabled: major radius Rc="
             << Corona.radius << " Rg, tube radius="
             << Corona.height << " Rg\n";
      }

      return true;
    }

    else
    {
      cout << "Error - '" << str << "' is not a key word!\n";
      cout << "Program aborted...\n\n\n";
      return false;
    }
  }
}


// return the number of detectors in the phi-direction
int Model::ScreenNumPhi()
{
  return SNumPhi;
}


// return the number of detectors in the theta-direction
int Model::ScreenNumTheta()
{
  return SNumTheta;
}


// return the number of spectroscopic channels
int Model::res()
{
  return resolut;
}


// Return the number of photons to sample
long int Model::PhotonNum()
{
  return PNum;
}


// return the number of scattering regions defined
int Model::GetScattNum()
{
  return EndScattList;
}


// return parameters of scattering region with index i
scatterer Model::GetCloud(int i)
{
  return Cloud[i];
}


// return the maximum expansion of all scattering regions
double Model::GetMaxExp()
{
  return MaxExp;
}


// return the maximum radial velocity of all scattering regions
double Model::GetVMax()
{
  return VMax;
}


// return the minimum photon wavelength defined
double Model::LambdaMin()
{
  return LamMin;
}


// return the maximum photon wavelength defined
double Model::LambdaMax()
{
  return LamMax;
}


// return the photon wavelength range covered
double Model::LambdaRange()
{
  return LamMax - LamMin;
}


int Model::GetRandomSeed()
{
  return RandomSeed;
}


bool Model::HasFixedRandomSeed()
{
  return RandomSeedFixed;
}


int Model::GetSpinCase()
{
  return SpinCase;
}


bool Model::HasSpinCase()
{
  return SpinCaseEnabled;
}


double Model::GetSpinISCO()
{
  return SpinISCO;
}


double Model::GetSkyAxisPA()
{
  return SkyAxisPA;
}

bool Model::HasInclination()
{
  return InclinationEnabled;
}

double Model::GetInclinationDeg()
{
  return InclinationDeg;
}

double Model::GetSourceRA()
{
  return SourceRA;
}

double Model::GetSourceDec()
{
  return SourceDec;
}

double Model::GetSourceDistanceKpc()
{
  return SourceDistanceKpc;
}

double Model::GetBlackHoleMassSolar()
{
  return BlackHoleMassSolar;
}


// return the X resolution defined
double Model::GetResX()
{
  return ResX;
}


// return the Y resolution defined
double Model::GetResY()
{
  return ResY;
}


// set the output file name to name
void Model::GetOutputFile(char name[300])
{
  if (strchr(outputfile,'/') || strchr(outputfile,'\\'))
  {
    strcpy(name,outputfile);
  }
  else
  {
    strcpy(name,"outputs/");
    strcat(name,outputfile);
    strcat(name,"/");
    strcat(name,outputfile);
  }
}


// return the minimum gab between scattering clouds
double Model::GetGap()
{
  return Gap;
}


// return 'true' if the model space is defined symmetrically with
// respect to the xy-plane and 'false' otherwise
bool Model::GetPlaneSymmetry()
{
  return planesym;
}


// return 'true' if the model space is defined symmetrically with
// respect to the z-axis and 'false' otherwise
bool Model::GetAxisSymmetry()
{
  return axissym;
}


// return the parameters of the continuum emission region
CTRegion Model::GetCTR()
{
  return CTR;
}




// return the parameters of the broad line emission region
BLRegion Model::GetBLR()
{
  return BLR;
}


// return the parameters of the narrow line emission region
NLRegion Model::GetNLR()
{
  return NLR;
}


// return 'true' if a corona region is defined in config.ini
bool Model::HasCorona()
{
  return UseCorona;
}


// return 'true' if photons should be emitted from the corona
bool Model::HasCoronaSource()
{
  return UseCorona && Corona.source;
}


// return true if photons should be emitted from the XRB disk seed source
bool Model::HasXRBSeedSource()
{
  return SeedDisk.source;
}


bool Model::ShouldSaveV()
{
  return SaveV;
}


bool Model::ShouldSaveComponents()
{
  return SaveComponents;
}


bool Model::ShouldSaveStokesFiles()
{
  return SaveStokesFiles;
}


bool Model::ShouldSaveMap()
{
  return SaveMap;
}


int Model::GetMapThetaIndex()
{
  return MapThetaIndex;
}


double Model::GetMapEnergyMin()
{
  return MapEnergyMin;
}


double Model::GetMapEnergyMax()
{
  return MapEnergyMax;
}


// return the parameters of the corona region
CoronaRegion Model::GetCorona()
{
  return Corona;
}


// return the parameters of the XRB disk seed photon source
XRBSeedDisk Model::GetXRBSeedDisk()
{
  return SeedDisk;
}
