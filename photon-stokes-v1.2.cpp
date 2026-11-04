// Implementation of all routines concerning a model photon. The
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


#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <stdio.h>    
#include <time.h>    
#include "photon-stokes-v1.2.h"                  
#include "corona-stokes-v1.2.h"
#include "disk-reflection-stokes-v1.2.h"
#include "disk-reflection-polarization-stokes-v1.2.h"
#include "randomc-stokes-v1.2.h"   

#if (defined (__BORLANDC__) || defined (_MSC_VER)) && ! defined(_WINDOWS_)
  #include <conio.h>                   // define getch() function
  #define _GETCH_DEFINED_
#endif               
                 
using namespace std;



  int seed = (int)time(0);            // random seed
  // choose one of the random number generators:
  CRandomMersenne RanGen(seed);       // make instance of random number generator



// set the photon parameters position, flying direction, Stokes vector, weight,
// and wavelength to the given values (*)
void Photon::SetPhotonData(double setpos[],double setdirect[][3],
                           double setI,double setQ,double setU,double setV,
                           double setW,double setlambda)
{
  int i;
  Equal(pos,setpos);
  for(i = 0; i <= 2; i++)
  {
    Equal(direct[i],setdirect[i]);
  }
  I = setI;
  Q = setQ;
  U = setU;
  V = setV;
  W = setW;
  lambda = setlambda;
  abs = false;
  XRBDiskHit = false;
  XRBDiskHitRadius = 0;
  XRBDiskReflectionNum = 0;
  XRBDiskHitNum = 0;
  CoronaScatterNum = 0;
  OriginFlag = 0;
  XRBDiskAbsorbed = false;
}


// reset photon to an initial state, located in a continuum or line
// emission region, being unpolarized, and flying toward a random
// direction
void Photon::Reset(Model M)
{

  // origin of the model space
  const double origin[3] = {0,0,0};

  // default frame of the model space, xyz coordinate system
  const double frame[3][3] = {{1,0,0},{0,1,0},{0,0,1}};

  int i,j,line;
  double R,phi,theta,min,max,r,mu,polarization;
  CTRegion CTR;
  BLRegion BLR;
  NLRegion NLR;

  // set the co-moving frame of the photon equal to the
  // default coordinate system of the model space
  for(i=0; i<=2; i++)
  {
    for(j=0; j<=2; j++)
    {
      direct[i][j] = frame[i][j];
    }
  }

  // sample phi and theta for the flight direction assuming isotropic emission
  phi = (RanGen.Random()) * 2*pi;

  do
  {
    theta = acos(1 - 2*(RanGen.Random()));
  }
  while ((theta == 0) || (theta == pi));


  // turn the photon's frame into the flight direction
  TurnDirection(2,phi);
  TurnDirection(1,theta-pi/2);

  if (M.HasXRBSeedSource())
  {
    OriginFlag = 1;
    ResetFromXRBSeedDisk(M);
  }
  else if (M.HasCoronaSource())
  {
    OriginFlag = 2;
    ResetFromCorona(M);
  }
  else
  {
    OriginFlag = 0;
    // decide if the photon to create is a continuum photon, a broad
    // line photon, or a narrow line photon
    CTR = M.GetCTR();
    BLR = M.GetBLR();
    NLR = M.GetNLR();
    r = RanGen.Random();
    if (r < NLR.rli)
    {
      line = 2;
    }
    else
    {
      if (r < NLR.rli+BLR.rli)
      {
        line = 1;
      }
      else
      {
        line = 0;
      }
    }
  
    switch (line)
    {

      // continuum photon, sample starting position inside the
      // cylinder CTR
      case 0:  
        if (CTR.R==0)
        {
          for(i=0; i<=2; i++)
          {
            pos[i] = origin[i];
          }
        }
        else
        {
          r = RanGen.Random();
          R = sqrt((CTR.R-CTR.b)*(CTR.R-CTR.b)+4*r*CTR.R*CTR.b);
          r = RanGen.Random();
          phi = 2*pi*r;
          pos[0] = R*cos(phi);
          pos[1] = R*sin(phi);
          r = RanGen.Random();
          pos[2] = CTR.a*(2*r-1);
        }

        // sample the photon wavelength according to an I ~ nu^(-alpha)
        // intensity spectrum; shift the wavelength if the emission
        // region is in motion
        r = RanGen.Random();
        if (CTR.alpha!=0)
        {
          min = pow(M.LambdaMin(),CTR.alpha);
          max = pow(M.LambdaMax(),CTR.alpha);
          lambda = pow(min+r*(max-min),1/CTR.alpha);
        }
        else
        {
          min = M.LambdaMin();
          max = M.LambdaMax();
          lambda = min*pow(max/min,r);
        }
        if (CTR.vr!=0 || CTR.vphi!=0 || CTR.vtheta!=0)
        {
          LambdaShift(false,CTR.vtype,CTR.vr,CTR.vphi,CTR.vtheta);
        }
        break;

      // broad line photon, sample starting position inside the
      // cylinder BLR and wavelength according to a lorentz profile;
      // shift the wavelength if the emission region is in motion
      case 1:  
        r = RanGen.Random();
        R = sqrt((BLR.R-BLR.b)*(BLR.R-BLR.b)+4*r*BLR.R*BLR.b);
        r = RanGen.Random();
        phi = 2*pi*r;
        pos[0] = R*cos(phi);
        pos[1] = R*sin(phi);
        r = RanGen.Random();
        pos[2] = BLR.a*(2*r-1);
        r = RanGen.Random();
        lambda = GetLineLambda(M,1,r);
        if (BLR.vr!=0 || BLR.vphi!=0 || BLR.vtheta!=0)
        {
          LambdaShift(false,BLR.vtype,BLR.vr,BLR.vphi,BLR.vtheta);
        }
        break;

      // narrow line photon, sample starting position inside the
      // double-cone NLR and wavelength according to a lorentz profile;
      // shift the wavelength if the emission region is in motion
      default:
        r = RanGen.Random();
        R = pow(r*pow(NLR.b,3)+(1-r)*pow(NLR.a,3),1.0/3);
        r = RanGen.Random();
        phi = 2*pi*r;
        r = RanGen.Random();
        theta = acos(1-r*(1-cos(NLR.theta)));
        CartCoord(pos,phi,theta);
        r = RanGen.Random();
        if (r<0.5)
        {
          SMult(pos,R,pos);
        }
        else
        {
          SMult(pos,-R,pos);
        }
        r = RanGen.Random();
        lambda = GetLineLambda(M,2,r);
        if (NLR.vr!=0 || NLR.vphi!=0 || NLR.vtheta!=0)
        {
          LambdaShift(false,NLR.vtype,NLR.vr,NLR.vphi,NLR.vtheta);
        }
        break;
    }
  }

  // initialize the Stokes vector to an unpolarized photon
  I = 1;
  V = 0;
  Q = 0;
  U = 0;

  if (M.HasXRBSeedSource() && M.GetXRBSeedDisk().polarization==1)
  {
    mu = direct[0][2];
    if (mu<0) mu = -mu;
    if (mu>1) mu = 1;
    polarization = 0.1171*(1-mu)/(1+3.582*mu);

    // direct[1] lies in the disk plane for disk-seed photons. In the
    // photon Stokes basis, negative Q therefore represents the electric
    // vector parallel to the projected disk plane.
    Q = -polarization;
  }

  // the weight factor may be used in the future, currently it
  // is set to unity for all photons
  W = 1;

  // the photon is outside any scattering region
  ScattFlag = ScattMax;

  // no scatterings occurred yet
  ScattNum = 0;

  // assume electron scattering cross-section
  CS = TCS;

  // the photon did not change scattering regions
  PathLength = 0;

  // the photon is not absorbed
  abs = false;

  // reset XRB disk-hit diagnostics
  XRBDiskHit = false;
  XRBDiskHitRadius = 0;
  XRBDiskReflectionNum = 0;
  XRBDiskHitNum = 0;
  CoronaScatterNum = 0;
  XRBDiskAbsorbed = false;
}


static double XRBSeedBlackbodyPDF(double energy,double tin)
{
  double x,expx;

  if (energy<=0 || tin<=0)
  {
    return 0;
  }

  x = energy/tin;
  if (x>80)
  {
    return 0;
  }
  if (x<1.0e-6)
  {
    return x;
  }

  expx = exp(x);
  return x*x/(expx-1);
}


static double SampleXRBSeedPowerLaw(XRBSeedDisk Disk,double r)
{
  double min,max,energy;

  if (Disk.emax<Disk.emin)
  {
    energy = Disk.emin;
    Disk.emin = Disk.emax;
    Disk.emax = energy;
  }
  if (Disk.emax==Disk.emin)
  {
    return Disk.emin;
  }

  if (Disk.gamma==1)
  {
    energy = Disk.emin*pow(Disk.emax/Disk.emin,r);
  }
  else
  {
    min = pow(Disk.emin,1-Disk.gamma);
    max = pow(Disk.emax,1-Disk.gamma);
    energy = pow(min+r*(max-min),1/(1-Disk.gamma));
  }
  return energy;
}


static double SampleXRBSeedBlackbody(XRBSeedDisk Disk)
{
  int i;
  double energy,pdf,max_pdf;

  if (Disk.emax<Disk.emin)
  {
    energy = Disk.emin;
    Disk.emin = Disk.emax;
    Disk.emax = energy;
  }
  if (Disk.emax==Disk.emin)
  {
    return Disk.emin;
  }

  // The maximum of x^2/(exp(x)-1) is about 0.648.
  // Use a slightly larger envelope for bounded rejection sampling.
  max_pdf = 0.7;
  for(i=0; i<10000; i++)
  {
    energy = Disk.emin+RanGen.Random()*(Disk.emax-Disk.emin);
    pdf = XRBSeedBlackbodyPDF(energy,Disk.tin);
    if (RanGen.Random()*max_pdf<=pdf)
    {
      return energy;
    }
  }

  return 0.5*(Disk.emin+Disk.emax);
}


static double SampleStandardNormal()
{
  double r1,r2;

  do
  {
    r1 = RanGen.Random();
  }
  while (r1<=0);
  r2 = RanGen.Random();

  return sqrt(-2*log(r1))*cos(2*pi*r2);
}


static void ApplyCoronaComptonShift(CoronaRegion Corona,double& lambda)
{
  double theta,energy,limit,meanFrac,sigmaLog,meanLog,amp,newEnergy;
  double recoil,gauss;

  if (!Corona.compton || Corona.kTe<=0 || lambda<=0)
  {
    return;
  }

  // Approximate thermal Comptonization for a single scattering.
  // The mean fractional gain follows the non-relativistic Kompaneets
  // estimate, with a recoil term and a stochastic thermal spread.
  energy = 12.3984193/lambda;
  if (energy<=0)
  {
    return;
  }

  theta = Corona.kTe/511.0;
  limit = 3*Corona.kTe;
  if (theta<=0 || limit<=0)
  {
    return;
  }

  recoil = energy/511.0;
  meanFrac = 4*theta+16*theta*theta-recoil;

  // Fade the up-scattering drift as photons approach the Wien limit.
  if (energy<limit && meanFrac>0)
  {
    meanFrac = meanFrac*(1-energy/limit);
  }

  sigmaLog = 0.5*sqrt(2*theta);
  if (sigmaLog>0.7)
  {
    sigmaLog = 0.7;
  }
  if (1+meanFrac<=0)
  {
    meanLog = -0.5;
  }
  else
  {
    meanLog = log(1+meanFrac)-0.5*sigmaLog*sigmaLog;
  }

  gauss = SampleStandardNormal();
  amp = exp(meanLog+sigmaLog*gauss);
  if (amp<0.05)
  {
    amp = 0.05;
  }
  if (amp>8.0)
  {
    amp = 8.0;
  }

  newEnergy = energy*amp;
  if (newEnergy>limit)
  {
    newEnergy = limit*(0.9+0.1*RanGen.Random());
  }
  if (newEnergy<=0)
  {
    return;
  }

  if (newEnergy!=energy)
  {
    lambda = 12.3984193/newEnergy;
  }
}


// reset photon to an initial state emitted by the XRB accretion disk
void Photon::ResetFromXRBSeedDisk(Model M)
{
  int i,j;
  double R,phi,mu,theta,energy;
  const double frame[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
  XRBSeedDisk Disk;

  Disk = M.GetXRBSeedDisk();

  if (Disk.rout<Disk.rin)
  {
    R = Disk.rin;
    Disk.rin = Disk.rout;
    Disk.rout = R;
  }
  if (Disk.rout<=0)
  {
    Disk.rout = 1;
  }
  if (Disk.rin<0)
  {
    Disk.rin = 0;
  }

  R = sqrt(Disk.rin*Disk.rin+
           RanGen.Random()*(Disk.rout*Disk.rout-Disk.rin*Disk.rin));
  phi = 2*pi*RanGen.Random();
  pos[0] = R*cos(phi);
  pos[1] = R*sin(phi);
  pos[2] = 0;

  // Emit into the upper hemisphere. The Lambert option samples
  // p(mu)=2*mu, where mu is the cosine to the disk normal.
  for(i=0; i<=2; i++)
  {
    for(j=0; j<=2; j++)
    {
      direct[i][j] = frame[i][j];
    }
  }
  phi = 2*pi*RanGen.Random();
  if (Disk.angular==1)
  {
    mu = sqrt(RanGen.Random());
  }
  else
  {
    mu = RanGen.Random();
  }
  theta = acos(mu);
  TurnDirection(2,phi);
  TurnDirection(1,theta-pi/2);

  if (Disk.emin<=0 || Disk.emax<=0)
  {
    lambda = 0.5*(M.LambdaMin()+M.LambdaMax());
    return;
  }
  if (Disk.emax<Disk.emin)
  {
    energy = Disk.emin;
    Disk.emin = Disk.emax;
    Disk.emax = energy;
  }
  if (Disk.spectrum==1)
  {
    energy = SampleXRBSeedBlackbody(Disk);
  }
  else
  {
    energy = SampleXRBSeedPowerLaw(Disk,RanGen.Random());
  }
  lambda = 12.3984193/energy;
}


// reset photon to an initial state emitted by the corona
void Photon::ResetFromCorona(Model M)
{
  double energy;
  CoronaRegion Corona;

  Corona = M.GetCorona();

  CoronaSamplePosition(Corona,RanGen.Random(),RanGen.Random(),
                       RanGen.Random(),pos);

  energy = GetCoronaEnergy(Corona,RanGen.Random());
  if (energy>0)
  {
    lambda = 12.3984193/energy;
  }
  else
  {
    lambda = 0.5*(M.LambdaMin()+M.LambdaMax());
  }
}


// sample a photon energy in keV from the corona power-law spectrum
double Photon::GetCoronaEnergy(CoronaRegion Corona,double r)
{
  return CoronaSampleEnergy(Corona,r);
}


// sample a wavelength for a line photon, decide according to the
// parameter 'line' if the photon is a broad or a narrow line
// photon. The wavelength is sampled assuming a Lorentz-profile of the
// intensity versus frequency spectrum
double Photon::GetLineLambda(Model M,int line,double r)
{
  double lambda,l1,l2,l,f;
  BLRegion BLR;
  NLRegion NLR;

  // broad line photon
  if (line == 1)
  {
    BLR = M.GetBLR();
    l1 = M.LambdaMin();
    l2 = M.LambdaMax();
    do
    {
      l = (l1+l2)/2;
      f = 0.5*log(pow((l-BLR.lambda),2)+pow(BLR.width,2))+
          BLR.lambda/BLR.width*atang(l-BLR.lambda,BLR.width)-
          BLR.low-r*BLR.norm;
      if (f >= 0)
        l2 = l;
      else
        l1 = l;
    }
    while (l2-l1>LinePrecision);
    lambda = (l1+l2)/2;
  }

  // narrow line photon
  else
  {
    NLR = M.GetNLR();
    l1 = M.LambdaMin();
    l2 = M.LambdaMax();
    do
    {
      l = (l1+l2)/2;
      f = 0.5*log(pow((l-NLR.lambda),2)+pow(NLR.width,2))+
          NLR.lambda/NLR.width*atang(l-NLR.lambda,NLR.width)-
          NLR.low-r*NLR.norm;
      if (f >= 0)
        l2 = l;
      else
        l1 = l;
    }
    while (l2-l1>LinePrecision);
    lambda = (l1+l2)/2;
  }

  return lambda;
}


// if the scattering region Cloud is in the photon's flight direction,
// then set 'tpre' to the distance to the cloud's closest border,
// otherwise leave 'tpre' unchanged
void Photon::FindNextBorder(scatterer Cloud,double& tpre,double pos[],
                            double direct[])
{
  double t;

  t = tpre;
  switch(Cloud.shape)
  {

    // check for a toroidal scattering region
    case 0:
      TorusCoinc(pos,direct,t,Cloud.a,Cloud.b,Cloud.R,Cloud.c);
      break;

    // check for a cylindrical scattering region, consisting of four
    // connected surfaces
    case 1:
      CylinderCoinc(pos,direct,t,Cloud.R-Cloud.b,Cloud.a,Cloud.c);
      CylinderCoinc(pos,direct,t,Cloud.R+Cloud.b,Cloud.a,Cloud.c);
      RingCoinc(pos,direct,t,Cloud.R-Cloud.b,Cloud.R+Cloud.b,-Cloud.a+Cloud.c);
      RingCoinc(pos,direct,t,Cloud.R-Cloud.b,Cloud.R+Cloud.b,Cloud.a+Cloud.c);
      break;

    // check for a double-conical scattering region, consisting of three
    // connected surfaces
    case 2:
      ConeCoinc(pos,direct,Cloud.a,Cloud.b,Cloud.theta,t);
      SphereSegCoinc(pos,direct,Cloud.pos,t,Cloud.a,Cloud.theta,true);
      SphereSegCoinc(pos,direct,Cloud.pos,t,Cloud.b,Cloud.theta,true);
      break;

    // check for a spherical scattering region
    case 3:
      SphereSegCoinc(pos,direct,Cloud.pos,t,Cloud.a,pi/2,true);	
      break;

    // check for a flared disc scattering region, consisting of three
    // connected surfaces
    default:
      ConeCoinc(pos,direct,Cloud.a,Cloud.b,pi/2-Cloud.theta,t);
      SphereSegCoinc(pos,direct,Cloud.pos,t,Cloud.a,Cloud.theta,false);
      SphereSegCoinc(pos,direct,Cloud.pos,t,Cloud.b,Cloud.theta,false);	
      break;
  }

  // find the closest scattering region
  if (t<tpre)
  {
    tpre = t;
  }

}


// if the corona is in the photon's flight direction, then set 'tpre'
// to the distance to the closest corona border, otherwise leave it unchanged
void Photon::FindNextCoronaBorder(CoronaRegion Corona,double& tpre,
                                  double pos[],double direct[])
{
  CoronaFindNextBorder(Corona,tpre,pos,direct);
}


// search for the nearest scattering region along the photons flight direction,
// if there is such cloud, then move the photon to the clouds border
bool Photon::ShiftToNextScatt(Model M)
{
  int i;
  double t,tpre,tcor,NewPos[3],OldPos[3];
  scatterer Cloud;
  CoronaRegion Corona;

  // set maximum flight distance that is possible
  tpre = M.GetMaxExp();
  t = M.GetMaxExp();
  Equal(OldPos,pos);

  // assume photon to be inside no scattering region
  ScattFlag = ScattMax;

  // loop over all defined scattering regions
  for(i=0; i<M.GetScattNum(); i++)
  {
    Cloud = M.GetCloud(i);

    // if the photon already is inside some scattering cloud
    if (InsideCloud(Cloud))
    {
      ScattFlag = i;

      // do not displace the photon
      t = -M.GetGap();

      break;
    }

    // search for the closest scattering region in the photons
    // flight direction
    else
    {
      FindNextBorder(Cloud,t,pos,direct[0]);

      // if the distance to the Cloud is smaller than the previous one
      // then re-define the distance
      if (t<tpre)
      {
        ScattFlag = i;
        tpre = t;
      }
    }
  }

  // include the optional corona as an electron scattering region
  if (M.HasCorona() && t>=0)
  {
    Corona = M.GetCorona();
    if (Corona.density>0 || Corona.tau>0)
    {
      if (InsideCorona(Corona,0))
      {
        ScattFlag = CoronaFlag;
        t = -M.GetGap();
      }
      else
      {
        tcor = tpre;
        FindNextCoronaBorder(Corona,tcor,pos,direct[0]);
        if (tcor<tpre)
        {
          ScattFlag = CoronaFlag;
          tpre = tcor;
          t = tcor;
        }
      }
    }
  }

  // new scattering region found
  if (ScattFlag < ScattMax)
  {
    double step;
    bool inside;

    // use the Thomson scattering cross-section (unit ccm/parsec)
    Cloud = M.GetCloud(ScattFlag);
    CS = TCS;

    // displace the photon, add 'Gap' to make sure that the photon is
    // declared to be inside the cloud
    step = M.GetGap();
    inside = false;
    for (i=0; i<12; i++)
    {
      SMult(NewPos,t+step,direct[0]);
      VecSum(NewPos,OldPos,NewPos);
      SetPosition(NewPos);
      if (InsideCloud(Cloud,M.GetGap()))
      {
        inside = true;
        break;
      }
      step = 0.1*step;
    }

    if (!inside)
    {
      SetPosition(OldPos);
      ScattFlag = ScattMax;
      return false;
    }

    return true;
  }

  // corona found
  else if (ScattFlag == CoronaFlag)
  {
    double step;
    bool inside;

    CS = TCS;

    step = M.GetGap();
    inside = false;
    for (i=0; i<12; i++)
    {
      SMult(NewPos,t+step,direct[0]);
      VecSum(NewPos,OldPos,NewPos);
      SetPosition(NewPos);
      if (InsideCorona(Corona,M.GetGap()))
      {
        inside = true;
        break;
      }
      step = 0.1*step;
    }

    if (!inside)
    {
      SetPosition(OldPos);
      ScattFlag = ScattMax;
      return false;
    }

    return true;
  }

  // no scattering region found
  else
  {
    return false;
  }

}


// sample a path length for the photon and shift it within the given
// scattering region. If the photon leaves the border of the region
// return false, otherwise return true
bool Photon::Proceed(Model M)
{
  double s,r,t,density,length,connect[3],OldPos[3],NewPos[3];
  scatterer Cloud;
  CoronaRegion Corona;

  // save the old photon position
  Equal(OldPos,pos);

  // sample a path length within the current scattering region
  if (ScattFlag == CoronaFlag)
  {
    Corona = M.GetCorona();
    density = Corona.density;
    if (density<=0 && Corona.tau>0)
    {
      if (Corona.type==0)
      {
        length = Corona.height;
      }
      else if (Corona.type==1 && Corona.height>0)
      {
        length = Corona.height;
      }
      else
      {
        length = Corona.radius;
      }
      if (length>0)
      {
        density = Corona.tau/(TCS*length);
      }
    }
  }
  else
  {
    Cloud = M.GetCloud(ScattFlag);
    density = Cloud.Density;
  }

  if (density<=0)
  {
    return false;
  }

  do
  {
    r = RanGen.Random();
  }
  while (r==0);

  if (PathLength > 0)
  {
    s = PathLength/(density*CS);
  }  

  // sample a path length
  else
  {
    do
    {
      r = double(rand())/RAND_MAX;
    }
    while (r==0);
    s = -log(r)/(density*CS);
  }

  // check if photon is leaving the current scattering region
  t = s;
  if (ScattFlag == CoronaFlag)
  {
    FindNextCoronaBorder(Corona,t,pos,direct[0]);
  }
  else
  {
    FindNextBorder(Cloud,t,pos,direct[0]);
  }

  // photon leaves the scattering region
  if (t < s)
  {

    // shift the photon just outside the current scattering region
    SMult(connect,t+M.GetGap(),direct[0]);
    VecSum(NewPos,pos,connect);
    if (ReflectFromXRBDisk(M,OldPos,NewPos))
    {
      PathLength = 0;
      return false;
    }
    Equal(pos,NewPos);

    // save remaining path length for next scattering region
    PathLength = (s-t)*density*CS;

    return false;
  }

  // photon stays in side the scattering region
  else
  {
    // displace the photon to the new position
    SMult(connect,s,direct[0]);
    VecSum(NewPos,pos,connect);
    if (ReflectFromXRBDisk(M,OldPos,NewPos))
    {
      PathLength = 0;
      return false;
    }
    Equal(pos,NewPos);

    // set path length to zero
    PathLength = 0.;

    return true;
  }

}


// return true if the photon is located inside the scattering region Cloud,
// return false otherwise;
bool Photon::InsideCloud(scatterer Cloud,double gap)
{
  double a,b,c,con[3],r,tn,eps,scale,rho2,z2,limit;
  bool inside;

  scale = Cloud.a;
  if (Cloud.b>scale) scale = Cloud.b;
  if (Cloud.R>scale) scale = Cloud.R;
  if (scale<=0) scale = 1;
  if (gap>0)
  {
    eps = gap;
    if (eps<1.0e-12*scale)
    {
      eps = 1.0e-12*scale;
    }
  }
  else
  {
    eps = 0;
  }

  // assume photon outside cloud
  inside = false;
  switch(Cloud.shape)
  {

    // check for a toroidal scattering region
    case 0:
      a = Cloud.a+eps;
      b = Cloud.b+eps;
      if ((pos[2]-Cloud.c)*(pos[2]-Cloud.c)/(a*a)+pow(Cloud.R-sqrt(pos[0]*pos[0]+pos[1]*pos[1]),2)/
         (b*b)<=1)
      {
        inside = true;
      }
      break;

    // check for a cylindrical scattering region
    case 1:
      r = sqrt(pos[0]*pos[0]+pos[1]*pos[1]);
      if ((r>=Cloud.R-Cloud.b-eps)&&(r<=Cloud.R+Cloud.b+eps)&&(fabs(pos[2]-Cloud.c)<=Cloud.a+eps))
      {
        inside = true;
      }
      break;

    // check for a double-conical scattering region
    case 2:
      r = SProd(pos,pos);
      tn = tan(pi/2-Cloud.theta);
      rho2 = pos[0]*pos[0]+pos[1]*pos[1];
      z2 = pos[2]*pos[2];
      limit = eps*(sqrt(r)+1);
      a = Cloud.a-eps;
      if (a<0) a = 0;
      b = Cloud.b+eps;
      if ((tn*tn*rho2<=z2+limit)&&
          (r>=a*a)&&(r<=b*b))
      {
        inside = true;
      }
      break;

    // check for a spherical scattering region
    case 3:
      VecDiff(con,pos,Cloud.pos);
      a = Cloud.a+eps;
      if (SProd(con,con) <= a*a)
      {
        inside = true;
      }
      break;

    // check for a flared disc scattering region
    default:
      r = SProd(pos,pos);
      tn = tan(Cloud.theta);
      rho2 = pos[0]*pos[0]+pos[1]*pos[1];
      z2 = pos[2]*pos[2];
      limit = eps*(sqrt(r)+1);
      a = Cloud.a-eps;
      if (a<0) a = 0;
      b = Cloud.b+eps;
      c = tn*tn*rho2;
      if ((c+limit>=z2)&&
          (r>=a*a)&&(r<=b*b))
      {
        inside = true;
      }
      break;

  }
  return inside;
}


// return true if the photon is located inside the corona region
bool Photon::InsideCorona(CoronaRegion Corona,double gap)
{
  return CoronaInside(Corona,pos,gap);
}


// return true if the current scattering region is the corona
bool Photon::InCorona()
{
  return ScattFlag == CoronaFlag;
}


// return true if the photon would hit the XRB disk plane before detection
bool Photon::WouldHitXRBDisk(Model M,double& radius)
{
  double t,x,y;
  XRBSeedDisk Disk;

  radius = 0;
  Disk = M.GetXRBSeedDisk();

  if (Disk.reflectRout<=Disk.reflectRin || Disk.reflectRout<=0)
  {
    return false;
  }

  // Count upper-side illumination of a flat disk in the z=0 plane.
  if (pos[2]<=0 || direct[0][2]>=0)
  {
    return false;
  }

  t = -pos[2]/direct[0][2];
  if (t<=0)
  {
    return false;
  }

  x = pos[0]+t*direct[0][0];
  y = pos[1]+t*direct[0][1];
  radius = sqrt(x*x+y*y);

  return (radius>=Disk.reflectRin && radius<=Disk.reflectRout);
}


bool Photon::HasHitXRBDisk(double& radius)
{
  radius = XRBDiskHitRadius;
  return XRBDiskHit;
}


int Photon::GetXRBDiskReflectionNum()
{
  return XRBDiskReflectionNum;
}


int Photon::GetCoronaScatterNum()
{
  return CoronaScatterNum;
}


bool Photon::IsReflectedComponent()
{
  return XRBDiskReflectionNum>0;
}


bool Photon::IsDirectComponent()
{
  return XRBDiskReflectionNum==0;
}


bool Photon::AbsorbedByXRBDisk()
{
  return XRBDiskAbsorbed;
}


bool Photon::ReflectEscapedXRBDisk(Model M)
{
  double t,OldPos[3],NewPos[3];

  if (pos[2]<=0 || direct[0][2]>=0)
  {
    return false;
  }

  t = -pos[2]/direct[0][2];
  if (t<=0)
  {
    return false;
  }

  Equal(OldPos,pos);
  NewPos[0] = pos[0]+t*direct[0][0];
  NewPos[1] = pos[1]+t*direct[0][1];
  NewPos[2] = -M.GetGap();

  return ReflectFromXRBDisk(M,OldPos,NewPos);
}


bool Photon::ReflectFromXRBDisk(Model M,double OldPos[],double NewPos[])
{
  double Incoming[3];
  DiskReflectionPolarization Polarization;
  DiskReflectionEvent Event;

  Equal(Incoming,direct[0]);
  DiskReflectionApply(M,OldPos,NewPos,Incoming,direct[1],direct[2],
                      I,Q,U,V,Event);
  if (!Event.hit)
  {
    return false;
  }

  XRBDiskHitNum++;
  if (!XRBDiskHit)
  {
    XRBDiskHit = true;
    XRBDiskHitRadius = Event.radius;
  }

  if (Event.absorbed)
  {
    Equal(pos,Event.pos);
    abs = true;
    XRBDiskAbsorbed = true;
    return true;
  }
  if (!Event.reflected)
  {
    return false;
  }

  Equal(pos,Event.pos);
  SetDirectionFrame(Event.direct);
  Equal(Polarization.axis,Event.polarAxis);
  Polarization.Q = Event.Q;
  Polarization.U = Event.U;
  Polarization.V = Event.V;
  DiskReflectionRotateToFrame(Polarization,direct[1],direct[2],Q,U);
  V = Event.V;
  XRBDiskReflectionNum++;
  return true;
}


void Photon::SetDirectionFrame(double NewDirect[])
{
  double norm;

  Equal(direct[0],NewDirect);
  norm = sqrt(SProd(direct[0],direct[0]));
  if (norm<=0)
  {
    direct[0][0] = 0;
    direct[0][1] = 0;
    direct[0][2] = 1;
  }
  else
  {
    direct[0][0] = direct[0][0]/norm;
    direct[0][1] = direct[0][1]/norm;
    direct[0][2] = direct[0][2]/norm;
  }

  if (fabs(direct[0][2])<0.999)
  {
    direct[1][0] = -direct[0][1];
    direct[1][1] = direct[0][0];
    direct[1][2] = 0;
  }
  else
  {
    direct[1][0] = 1;
    direct[1][1] = 0;
    direct[1][2] = 0;
  }
  norm = sqrt(SProd(direct[1],direct[1]));
  direct[1][0] = direct[1][0]/norm;
  direct[1][1] = direct[1][1]/norm;
  direct[1][2] = direct[1][2]/norm;

  direct[2][0] = direct[0][1]*direct[1][2]-direct[0][2]*direct[1][1];
  direct[2][1] = direct[0][2]*direct[1][0]-direct[0][0]*direct[1][2];
  direct[2][2] = direct[0][0]*direct[1][1]-direct[0][1]*direct[1][0];
}


// turn the photon's flight direction by the angle 'alpha' and around the
// coordinate axis given by int; mathematically this is a matrix multiplication
// of the rotation matrix mat with each unity vector direct[i] of the photon's
// co-moving coordinate system
void Photon::TurnDirection(int axis,double ep)
{
  int i,j;
  double sn,cs,temp,norm;
  double tempdirect[3];
  double mat[3][3];

  sn = sin(ep);
  cs = cos(ep);

  // prepare the rotation matrix mat
  for(i = 0; i <=2; i++) mat[i][i] = pow(direct[axis][i],2) * (1 - cs) + cs;
  temp = direct[axis][0] * direct[axis][1] * (1 - cs);
  mat[0][1] = temp - direct[axis][2] * sn;
  mat[1][0] = temp + direct[axis][2] * sn;
  temp = direct[axis][0] * direct[axis][2] * (1 - cs);
  mat[0][2] = temp + direct[axis][1] * sn;
  mat[2][0] = temp - direct[axis][1] * sn;
  temp = direct[axis][1] * direct[axis][2] * (1 - cs);
  mat[1][2] = temp - direct[axis][0] * sn;
  mat[2][1] = temp + direct[axis][0] * sn;

  // matrix multiplication and re-normalization of the unity vectors
  // to avoid numerical inaccuracies
  for(i = 0; i <= 2; i++)
  {
    MatVecMult(mat,direct[i],tempdirect);
    for(j = 0; j <= 2; j++) direct[i][j] = tempdirect[j];
    norm = sqrt(SProd(direct[i],direct[i]));
    for(j = 0; j <= 2; j++) direct[i][j] = 1/norm * direct[i][j];
  }

}


// turn the photon's polar plane by the angle 'alpha' around the axis given by
// its flight direction; this is a rotation of the polarization plane, defined
// by direct[1] and direct[2], within itself; the Stokes vector is adjusted to
// keep the correct polarization
void Photon::TurnPolarPlane(double alpha)
{
  double QNew,UNew,c1,c2;

  c1 = cos(2*alpha);
  c2 = sin(2*alpha);
  QNew = c1*Q - c2*U;
  UNew = c2*Q + c1*U;
  Q = QNew;
  U = UNew;
  TurnDirection(0,alpha);
}


// perform an electron scattering event and change all affected
// photon parameters
void Photon::ElectronScattering(Model M)
{
  double PL,theta,phi;
  CoronaRegion Corona;
  scatterer Cloud;

  if (ScattFlag == CoronaFlag)
  {
    Corona = M.GetCorona();
    Cloud.vtype = 0;
    Cloud.vr = 0;
    Cloud.vphi = 0;
    Cloud.vtheta = 0;
  }
  else
  {
    Cloud = M.GetCloud(ScattFlag);
  }

  // if the scattering cloud is in motion transform
  // the photon wavelength
  if (Cloud.vr!=0 || Cloud.vphi!=0 || Cloud.vtheta!=0)
  {
    LambdaShift(true,Cloud.vtype,Cloud.vr,Cloud.vphi,Cloud.vtheta);
  }

  // sample a scattering angle theta
  theta = GetElecScattTheta(RanGen.Random());

  // compute the incident polarization and sample a scattering angle phi
  PL = sqrt(Q*Q+U*U)/I;
  phi = GetElecScattPhi(RanGen.Random(),PL,theta);

  // adjust the photon's flight direction and Stokes parameters
  TurnPolarPlane(pi-0.5*atang(U,Q)+phi);
  ElectTurnTheta(theta);

  // if the scattering cloud is in motion transform
  // the photon wavelength
  if (Cloud.vr!=0 || Cloud.vphi!=0 || Cloud.vtheta!=0)
  {
    LambdaShift(false,Cloud.vtype,Cloud.vr,Cloud.vphi,Cloud.vtheta);
  }
  if (ScattFlag == CoronaFlag)
  {
    ApplyCoronaComptonShift(Corona,lambda);
    CoronaScatterNum++;
  }

  // increase the number of scatterings
  ScattNum = ScattNum + 1;

}


// sample an electron scattering angle theta according to the
// random number r; the necessary trigonometric equation is
// solved by an iterative method
double  Photon::GetElecScattTheta(double r)
{
  double p1,p2,p,f;

  p1 = 0;
  p2 = pi;
  do
  {
    p = (p2+p1)/2;
    f = 0.5+0.125*cos(p)*(pow(sin(p),2)-4) - r;
    if (f >= 0)
      p2 = p;
    else
      p1 = p;
  }
  while (p2-p1 > pi/180*PFPrecision);
  return (p2+p1)/2;
}


// sample an electron scattering angle phi according to the
// random number r, incident polarization PL, and the scattering
// angle theta; the necessary trigonometric equation is solved by
// an iterative method
double Photon::GetElecScattPhi(double r,double PL,double theta)
{
  double p1,p2,p,f,N,csq;

  p1 = 0;
  p2 = 2*pi;
  N = 1/(2*pi);
  csq = pow(cos(theta),2);
  do
  {
    p = (p2+p1)/2;
    f = N*(p - (1-csq)/(1+csq)*PL*sin(2*p)/2) - r;
    if (f >= 0)
      p2 = p;
    else
      p1 = p;
  }
  while (p2-p1 > pi/180*PFPrecision);
  return (p2+p1)/2;
}


// turn the flight direction of the photon due to electron scattering by
// the angle theta, adjust the Stokes parameters
void Photon::ElectTurnTheta(double theta)
{
  double INew,VNew,QNew,UNew,cs;

  cs = cos(theta);


  // compute new Stokes parameters an normalize them to I=1
  INew = 0.5*(1+cs*cs)*I + 0.5*(cs*cs-1)*Q;
  QNew = 0.5*(cs*cs-1)*I + 0.5*(1+cs*cs)*Q;
  UNew = cs * U;
  VNew = cs * V;
  if (INew!=0)
  {
    Q = QNew/INew;
    U = UNew/INew;
    V = VNew/INew;
    I = 1;
  }

  // turn the flight direction by the scattering angle theta
  TurnDirection(2,theta);

}


// Transform the photon wavelength to (ToFrame=true) or from (ToFrame=false)
// a reference frame moving in the direction of the velocity (vr,vphi,vtheta);
// the new wavelength is computed by a relativistic Doppler shift formula.
void Photon::LambdaShift(bool ToFrame,int vtype,double vr,double vphi,double vtheta)
{
  double betak,beta2;
  double PhiP, ThetaP;
  double v[3],p[3],vroh,vz;

  // by decommenting and adjusting these lines, radial dependecies of 
  // the velocity can be assigned
  // const double r0 = 1.;
  // const double vr0 = 1.;
  // const double vphi0 = 1.;
  // const double vtheta0 = 1.;
  // const double alpha_r = -2;
  // const double alpha_p = -2;
  // const double alpha_t = -2;
  // vtype = 0;
  // vr = vr0 * pow(GetDistance()/r0,alpha_r);
  // vphi = vphi0 * pow(GetDistance()/r0,alpha_p);
  // vtheta = vtheta0 * pow(GetDistance()/r0,alpha_t);

  // transform the velocity from local to global Cartesian coordinates
  SMult(p,1/GetDistance(),pos);
  SphericalCoords(p,PhiP,ThetaP); 
  switch(vtype)
  {
    
    //input velocity given in cylindrical coordinates
    case 1: 
      vroh = vr;
      vz   = vtheta;
      break;

    //input velocity given in spherical coordinates
    default: 
      vroh = vr*sin(ThetaP)-vtheta*cos(ThetaP);
      vz   = vr*cos(ThetaP)+vtheta*sin(ThetaP);
      break;

  }  
  v[0] = vroh*cos(PhiP) - vphi*sin(PhiP);
  v[1] = vroh*sin(PhiP) + vphi*cos(PhiP);
  v[2] = vz; 

  // project onto the photon velocity and compute beta^2
  betak = SProd(v,direct[0])/vc;
  beta2 = SProd(v,v)/(vc*vc);

  // compute the wavelength shift
  if (ToFrame)
  {
    lambda = lambda*(1+betak)/sqrt(1-beta2);
  }
  else
  {
    lambda = lambda*(1-betak)/sqrt(1-beta2);
  }

}


// adjust the Stokes parameters of the photon to the polarization plane of
// a detector in the direction given by theta and phi
void Photon::AdjustFrame(double phi,double theta)
{
  double vp[3],sp,angle,FDirect[3];

  // find the unity vector FDirect of the detector's polarization plane
  // that points 'down-wards'
  if (theta <= pi/2)
  {
    CartCoord(FDirect,phi,theta+pi/2);
  }
  else
  {
    CartCoord(FDirect,phi+pi,3*pi/2-theta);
  }

  // compute the angle between FDirect and the co-moving unity vector
  // direct[2]; from their scalar product sp; check if sp is within the
  // allowed limits
  sp = SProd(FDirect,direct[2]);
  if (sp > 1)
  {
    sp = 1;
  }
  if (sp < -1)
  {
    sp = -1;
  }
  angle = acos(sp);

  // compute the vector product between FDirect and direct[2] to
  // decide in which direction the photon's polarization plane has
  // to be turned
  vp[0] = FDirect[1]*direct[2][2]-direct[2][1]*FDirect[2];
  vp[1] = FDirect[2]*direct[2][0]-direct[2][2]*FDirect[0];
  vp[2] = FDirect[0]*direct[2][1]-direct[2][0]*FDirect[1];
  if (SProd(vp,direct[0]) >= 0)
  {
    TurnPolarPlane(-angle);
  }
  else
  {
    TurnPolarPlane(angle);
  }
}


// return the photon wavelength
double Photon::GetLambda()
{
  return lambda;
}


// set photon position to 'setpos'
void Photon::SetPosition(double setpos[])
{
  int i;
  for(i = 0; i <= 2; i++)
    pos[i] = setpos[i];
}


// set 'getpos' to the current photon position
void Photon::GetPosition(double getpos[])
{
  int i;
  for(i = 0; i <= 2; i++)
    getpos[i] = pos[i];
}


// return the x-component of the flight direction
double Photon::GetDirect1()
{
  return direct[0][0];
}


// return the y-component of the flight direction
double Photon::GetDirect2()
{
  return direct[0][1];
}


// return the z-component of the flight direction
double Photon::GetDirect3()
{
  return direct[0][2];
}


// return the x-component of the photon position
double Photon::GetPos1()
{
  return pos[0];
}


// return the y-component of the photon position
double Photon::GetPos2()
{
  return pos[1];
}


// return the z-component of the photon position
double Photon::GetPos3()
{
  return pos[2];
}


// set the scattering region number to 'NewFlag'
void Photon::SetScattFlag(int NewFlag)
{
  ScattFlag = NewFlag;
}


// return the current scattering region number
int Photon::GetFlag()
{
  return ScattFlag;
}


// return the number of scatterings the photon has undergone (*)
int Photon::GetScattNum()
{
  return ScattNum;
}


// return the distance of the photon from the origin of the model space
double Photon::GetDistance()
{
  return sqrt(SProd(pos,pos));
}


// return the I value of the Stokes vector
double Photon::GetI()
{
  return I;
}


// return the Q value of the Stokes vector
double Photon::GetQ()
{
  return Q;
}


// return the U value of the Stokes vector
double Photon::GetU()
{
  return U;
}


// return the V value of the Stokes vector
double Photon::GetV()
{
  return V;
}


// return the weight factor of the photon (*)
double Photon::GetW()
{
  return W;
}


// return the value of the photon's absorption flag
bool Photon::Absorbed()
{
  return abs;
}


// set the value of the photon's absorption flag to true
void Photon::SetAbsorbed()
{
  abs = true;
}










// Mesrenne Twister functions
void CRandomMersenne::Init0(int seed) {
   // Seed generator
   const uint32_t factor = 1812433253UL;
   mt[0]= seed;
   for (mti=1; mti < MERS_N; mti++) {
      mt[mti] = (factor * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
   }
}

void CRandomMersenne::RandomInit(int seed) {
   // Initialize and seed
   Init0(seed);

   // Randomize some more
   for (int i = 0; i < 37; i++) BRandom();
}


void CRandomMersenne::RandomInitByArray(int const seeds[], int NumSeeds) {
   // Seed by more than 32 bits
   int i, j, k;

   // Initialize
   Init0(19650218);

   if (NumSeeds <= 0) return;

   // Randomize mt[] using whole seeds[] array
   i = 1;  j = 0;
   k = (MERS_N > NumSeeds ? MERS_N : NumSeeds);
   for (; k; k--) {
      mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL)) + (uint32_t)seeds[j] + j;
      i++; j++;
      if (i >= MERS_N) {mt[0] = mt[MERS_N-1]; i=1;}
      if (j >= NumSeeds) j=0;}
   for (k = MERS_N-1; k; k--) {
      mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL)) - i;
      if (++i >= MERS_N) {mt[0] = mt[MERS_N-1]; i=1;}}
   mt[0] = 0x80000000UL;  // MSB is 1; assuring non-zero initial array

   // Randomize some more
   mti = 0;
   for (int i = 0; i <= MERS_N; i++) BRandom();
}


uint32_t CRandomMersenne::BRandom() {
   // Generate 32 random bits
   uint32_t y;

   if (mti >= MERS_N) {
      // Generate MERS_N words at one time
      const uint32_t LOWER_MASK = (1LU << MERS_R) - 1;       // Lower MERS_R bits
      const uint32_t UPPER_MASK = 0xFFFFFFFF << MERS_R;      // Upper (32 - MERS_R) bits
      static const uint32_t mag01[2] = {0, MERS_A};

      int kk;
      for (kk=0; kk < MERS_N-MERS_M; kk++) {    
         y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
         mt[kk] = mt[kk+MERS_M] ^ (y >> 1) ^ mag01[y & 1];}

      for (; kk < MERS_N-1; kk++) {    
         y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
         mt[kk] = mt[kk+(MERS_M-MERS_N)] ^ (y >> 1) ^ mag01[y & 1];}      

      y = (mt[MERS_N-1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
      mt[MERS_N-1] = mt[MERS_M-1] ^ (y >> 1) ^ mag01[y & 1];
      mti = 0;
   }
   y = mt[mti++];

   // Tempering (May be omitted):
   y ^=  y >> MERS_U;
   y ^= (y << MERS_S) & MERS_B;
   y ^= (y << MERS_T) & MERS_C;
   y ^=  y >> MERS_L;

   return y;
}


double CRandomMersenne::Random() {
   // Output random float number in the interval 0 <= x < 1
   // Multiply by 2^(-32)
   return (double)BRandom() * (1./(65536.*65536.));
}


int CRandomMersenne::IRandom(int min, int max) {
   // Output random integer in the interval min <= x <= max
   // Relative error on frequencies < 2^-32
   if (max <= min) {
      if (max == min) return min; else return 0x80000000;
   }
   // Multiply interval with random and truncate
   int r = int((double)(uint32_t)(max - min + 1) * Random() + min); 
   if (r > max) r = max;
   return r;
}


int CRandomMersenne::IRandomX(int min, int max) {
   // Output random integer in the interval min <= x <= max
   // Each output value has exactly the same probability.
   // This is obtained by rejecting certain bit values so that the number
   // of possible bit values is divisible by the interval length
   if (max <= min) {
      if (max == min) return min; else return 0x80000000;
   }
#ifdef  INT64_SUPPORTED
   // 64 bit integers available. Use multiply and shift method
   uint32_t interval;                    // Length of interval
   uint64_t longran;                     // Random bits * interval
   uint32_t iran;                        // Longran / 2^32
   uint32_t remainder;                   // Longran % 2^32

   interval = uint32_t(max - min + 1);
   if (interval != LastInterval) {
      // Interval length has changed. Must calculate rejection limit
      // Reject when remainder >= 2^32 / interval * interval
      // RLimit will be 0 if interval is a power of 2. No rejection then
      RLimit = uint32_t(((uint64_t)1 << 32) / interval) * interval - 1;
      LastInterval = interval;
   }
   do { // Rejection loop
      longran  = (uint64_t)BRandom() * interval;
      iran = (uint32_t)(longran >> 32);
      remainder = (uint32_t)longran;
   } while (remainder > RLimit);
   // Convert back to signed and return result
   return (int32_t)iran + min;

#else
   // 64 bit integers not available. Use modulo method
   uint32_t interval;                    // Length of interval
   uint32_t bran;                        // Random bits
   uint32_t iran;                        // bran / interval
   uint32_t remainder;                   // bran % interval

   interval = uint32_t(max - min + 1);
   if (interval != LastInterval) {
      // Interval length has changed. Must calculate rejection limit
      // Reject when iran = 2^32 / interval
      // We can't make 2^32 so we use 2^32-1 and correct afterwards
      RLimit = (uint32_t)0xFFFFFFFF / interval;
      if ((uint32_t)0xFFFFFFFF % interval == interval - 1) RLimit++;
   }
   do { // Rejection loop
      bran = BRandom();
      iran = bran / interval;
      remainder = bran % interval;
   } while (iran >= RLimit);
   // Convert back to signed and return result
   return (int32_t)remainder + min;

#endif
}









/***********************************************************************
                     End of program
***********************************************************************/

void EndOfProgram() {
  // This function takes care of whatever is necessary to do when the 
  // program is finished

  // It may be necessary to wait for the user to press a key
  // in order to prevent the output window from disappearing.
  // Remove the #ifdef and #endif lines to unconditionally wait for a key press;
  // Remove all three lines to not wait:
  #ifdef _GETCH_DEFINED_
  _getch();                            // wait for user to press a key
  #endif

  // It may be necessary to end the program with a linefeed:
  #if defined (__unix__) || defined (_MSC_VER)
  printf("\n");                        // end program with a linefeed
  #endif
}


/***********************************************************************
                     Error message function
***********************************************************************/

void FatalError(const char * ErrorText) {
  // This function outputs an error message and aborts the program.

  // Important: There is no universally portable way of outputting an 
  // error message. You may have to modify this function to output
  // the error message in a way that is appropriate for your system.


  // Check if FatalAppExit exists (this macro is defined in winbase.h)
  #ifdef FatalAppExit  

  // in Windows, use FatalAppExit:
  FatalAppExit(0, ErrorText);

  #else

  // in console mode, print error message  
  printf ("\n%s\n", ErrorText);
  EndOfProgram();

  #endif

  // Terminate program with error code
  exit(1);
}


