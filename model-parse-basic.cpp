#include "model-stokes-v1.2.h"
#include <istream>
#include <cstdlib>
#include <cstring>
#include <cstdio>

using namespace std;

bool Model::ParseGeneralOption(std::istream& in_stream,char str[200])
{
  bool accept = false;
  if (strcmp(str,"PhotonNum")==0)
  {
    in_stream >> str;
    PNum = atol(str);
    accept = true;
  }
  if (strcmp(str,"RandomSeed")==0)
  {
    in_stream >> str;
    RandomSeed = atoi(str);
    RandomSeedFixed = (RandomSeed>0);
    accept = true;
  }
  if (strcmp(str,"SpinCase")==0)
  {
    in_stream >> str;
    SpinCase = atoi(str);
    if (SpinCase==0)
    {
      SpinISCO = 6.0;
      SpinCaseEnabled = true;
      accept = true;
    }
    else if (SpinCase==1)
    {
      SpinISCO = 1.24;
      SpinCaseEnabled = true;
      accept = true;
    }
    else
    {
      accept = false;
    }
  }
  if (strcmp(str,"SkyAxisPA")==0)
  {
    in_stream >> str;
    SkyAxisPA = atof(str);
    while (SkyAxisPA < -180.0) SkyAxisPA += 360.0;
    while (SkyAxisPA >= 180.0) SkyAxisPA -= 360.0;
    accept = true;
  }
  if (strcmp(str,"InclinationDeg")==0)
  {
    in_stream >> str;
    InclinationDeg = atof(str);
    InclinationEnabled = true;
    accept = true;
  }
  if (strcmp(str,"SourceRA")==0)
  {
    in_stream >> str;
    SourceRA = atof(str);
    accept = true;
  }
  if (strcmp(str,"SourceDec")==0)
  {
    in_stream >> str;
    SourceDec = atof(str);
    accept = true;
  }
  if (strcmp(str,"SourceDistanceKpc")==0)
  {
    in_stream >> str;
    SourceDistanceKpc = atof(str);
    accept = true;
  }
  if (strcmp(str,"BlackHoleMassSolar")==0)
  {
    in_stream >> str;
    BlackHoleMassSolar = atof(str);
    accept = true;
  }
  if (strcmp(str,"LambdaMin")==0)
  {
    in_stream >> str;
    LamMin = atof(str);
    accept = true;
  }
  if (strcmp(str,"LambdaMax")==0)
  {
    in_stream >> str;
    LamMax = atof(str);
    accept = true;
  }
  return accept;
}

bool Model::ParseLegacyOption(std::istream& in_stream,char str[200])
{
  bool accept = false;
  int i,j;
  if (strcmp(str,"ContSource")==0)
  {
    in_stream >> str;
    CTR.R = atof(str);
    in_stream >> str;
    CTR.a = atof(str);
    in_stream >> str;
    CTR.b = atof(str);
    in_stream >> str;
    CTR.c = atof(str);
    in_stream >> str;
    CTR.alpha = atof(str);
	  in_stream >> str;
	  CTR.vtype = atoi(str);
	  in_stream >> str;
	  CTR.vr = atof(str);
	  in_stream >> str;
	  CTR.vphi = atof(str);
	  in_stream >> str;
	  CTR.vtheta = atof(str);
    accept = true;
  }
  if (strcmp(str,"BLSource")==0)
  {
    in_stream >> str;
    BLR.R = atof(str);
    in_stream >> str;
    BLR.a = atof(str);
    in_stream >> str;
    BLR.b = atof(str);
    in_stream >> str;
    BLR.lambda = atof(str);
    in_stream >> str;
    BLR.width = 0.5*atof(str);
    in_stream >> str;
    BLR.rli = 0.01*atof(str);
	  in_stream >> str;
	  BLR.vtype = atoi(str);
	  in_stream >> str;
	  BLR.vr = atof(str);
	  in_stream >> str;
	  BLR.vphi = atof(str);
	  in_stream >> str;
	  BLR.vtheta = atof(str);
    accept = true;
  }
  if (strcmp(str,"NLSource")==0)
  {
    in_stream >> str;
    NLR.a = atof(str);
    in_stream >> str;
    NLR.b = atof(str);
    in_stream >> str;
    NLR.theta = pi/180*atof(str);
    in_stream >> str;
    NLR.lambda = atof(str);
    in_stream >> str;
    NLR.width = 0.5*atof(str);
    in_stream >> str;
    NLR.rli = 0.01*atof(str);
	  in_stream >> str;
	  NLR.vtype = atoi(str);
	  in_stream >> str;
	  NLR.vr = atof(str);
	  in_stream >> str;
	  NLR.vphi = atof(str);
	  in_stream >> str;
	  NLR.vtheta = atof(str);
    accept = true;
  }
  if (strcmp(str,"Torus")==0)
  {
    i = EndScattList;
    Cloud[i].shape = 0;
    in_stream >> str;
    Cloud[i].R = atof(str);
    in_stream >> str;
    Cloud[i].a = atof(str);
    in_stream >> str;
    Cloud[i].b = atof(str);
    in_stream >> str;
    Cloud[i].c = atof(str);
    in_stream >> str;
    Cloud[i].Density = atof(str);
	  in_stream >> str;
	  Cloud[i].vtype = atoi(str);
    in_stream >> str;
    Cloud[i].vr = atof(str);
    in_stream >> str;
    Cloud[i].vphi = atof(str);
    in_stream >> str;
    Cloud[i].vtheta = atof(str);
    EndScattList++;
    accept = true;
  }
  if (strcmp(str,"Cylinder")==0)
  {
    i = EndScattList;
    Cloud[i].shape = 1;
    in_stream >> str;
    Cloud[i].R = atof(str);
    in_stream >> str;
    Cloud[i].a = atof(str);
    in_stream >> str;
    Cloud[i].b = atof(str);
    in_stream >> str;
    Cloud[i].c = atof(str);
    in_stream >> str;
    Cloud[i].Density = atof(str);
	  in_stream >> str;
	  Cloud[i].vtype = atoi(str);
    in_stream >> str;
    Cloud[i].vr = atof(str);
    in_stream >> str;
    Cloud[i].vphi = atof(str);
    in_stream >> str;
    Cloud[i].vtheta = atof(str);
    EndScattList++;
    accept = true;
  }
  if (strcmp(str,"Double-cone")==0)
  {
    i = EndScattList;
    Cloud[i].shape = 2;
    in_stream >> str;
    Cloud[i].a = atof(str);
    in_stream >> str;
    Cloud[i].b = atof(str);
    in_stream >> str;
    Cloud[i].theta = pi/180*atof(str);
    in_stream >> str;
    Cloud[i].Density = atof(str);
	  in_stream >> str;
	  Cloud[i].vtype = atoi(str);
    in_stream >> str;
    Cloud[i].vr = atof(str);
    in_stream >> str;
    Cloud[i].vphi = atof(str);
    in_stream >> str;
    Cloud[i].vtheta = atof(str);
    EndScattList++;
    accept = true;
  }
  if (strcmp(str,"Sphere")==0)
  {
    i = EndScattList;
    Cloud[i].shape = 3;
    for (j=0; j<=2; j++)
    {
      in_stream >> str;
      Cloud[i].pos[j] = atof(str);
    }
    in_stream >> str;
    Cloud[i].a = atof(str);
    in_stream >> str;
    Cloud[i].Density = atof(str);
	  in_stream >> str;
	  Cloud[i].vtype = atoi(str);
    in_stream >> str;
    Cloud[i].vr = atof(str);
    in_stream >> str;
    Cloud[i].vphi = atof(str);
    in_stream >> str;
    Cloud[i].vtheta = atof(str);
    EndScattList++;
    accept = true;
  }
  if (strcmp(str,"Flared-disk")==0)
  {
    i = EndScattList;
    Cloud[i].shape = 4;
    in_stream >> str;
    Cloud[i].a = atof(str);
    in_stream >> str;
    Cloud[i].b = atof(str);			
    in_stream >> str;
    Cloud[i].theta = pi/180*atof(str);
    in_stream >> str;
    Cloud[i].Density = atof(str);
	  in_stream >> str;
	  Cloud[i].vtype = atoi(str);
    in_stream >> str;
    Cloud[i].vr = atof(str);
    in_stream >> str;
    Cloud[i].vphi = atof(str);
    in_stream >> str;
    Cloud[i].vtheta = atof(str);
    EndScattList++;
    accept = true;
  }
  return accept;
}

bool Model::ParseOutputOption(std::istream& in_stream,char str[200])
{
  bool accept = false;
  int FileNum;
  char MStr[4];
  if (strcmp(str,"Gap")==0)
  {
    in_stream >> str;
    Gap = atof(str);
    accept = true;
  }
  if (strcmp(str,"PhiViewAng")==0)
  {
    in_stream >> str;
    SNumPhi = atoi(str);
    accept = true;
  }
  if (strcmp(str,"ThetaViewAng")==0)
  {
    in_stream >> str;
    SNumTheta = atoi(str);
    accept = true;
  }
  if (strcmp(str,"PlaneSym")==0)
  {
    in_stream >> str;
    planesym = (strcmp(str,"yes")==0);
    accept = true;
  }
  if (strcmp(str,"AxisSym")==0)
  {
    in_stream >> str;
    axissym = (strcmp(str,"yes")==0);
    accept = true;
  }
  if (strcmp(str,"SpectRes")==0)
  {
    in_stream >> str;
    resolut = atoi(str);
    accept = true;
  }
  if (strcmp(str,"ResX")==0)
  {
    in_stream >> str;
    ResX = atoi(str);
    accept = true;
  }
  if (strcmp(str,"ResY")==0)
  {
    in_stream >> str;
    ResY = atoi(str);
    accept = true;
  }
  if (strcmp(str,"SaveV")==0)
  {
    in_stream >> str;
    SaveV = (strcmp(str,"yes")==0);
    if (SaveV || strcmp(str,"no")==0)
    {
      accept = true;
    }
  }
  if (strcmp(str,"SaveComponents")==0)
  {
    in_stream >> str;
    SaveComponents = (strcmp(str,"yes")==0);
    if (SaveComponents || strcmp(str,"no")==0)
    {
      accept = true;
    }
  }
  if (strcmp(str,"SaveStokesFiles")==0)
  {
    in_stream >> str;
    SaveStokesFiles = (strcmp(str,"yes")==0);
    if (SaveStokesFiles || strcmp(str,"no")==0)
    {
      accept = true;
    }
  }
  if (strcmp(str,"SaveMap")==0)
  {
    in_stream >> str;
    SaveMap = (strcmp(str,"yes")==0);
    if (SaveMap || strcmp(str,"no")==0)
    {
      accept = true;
    }
  }
  if (strcmp(str,"MapThetaIndex")==0)
  {
    in_stream >> str;
    MapThetaIndex = atoi(str);
    accept = true;
  }
  if (strcmp(str,"MapEnergyRange")==0)
  {
    in_stream >> str;
    MapEnergyMin = atof(str);
    in_stream >> str;
    MapEnergyMax = atof(str);
    accept = true;
  }
	if (strcmp(str,"OutputFile")==0)
	{
	  in_stream >> outputfile;
	  in_stream >> str;
	  FileNum = atoi(str);
	  sprintf(MStr,"%d",FileNum);
	  strcat(outputfile,MStr);
	  accept = true;
	}
  return accept;
}
