#include "model-stokes-v1.2.h"
#include <istream>
#include <cstdlib>
#include <cstring>

using namespace std;

bool Model::ParseCoronaOption(std::istream& in_stream,char str[200],
                              bool& ReadCoronaSource)
{
  bool accept = false;
  int j;
  if (strcmp(str,"CoronaType")==0)
  {
    in_stream >> str;
    accept = true;
    if (strcmp(str,"slab")==0)
    {
      Corona.type = 0;
    }
    else if (strcmp(str,"toroidal")==0 || strcmp(str,"torus")==0 ||
             strcmp(str,"sphere")==0)
    {
      Corona.type = 1;
    }
    else if (strcmp(str,"lamp-post")==0)
    {
      Corona.type = 2;
    }
    else
    {
      accept = false;
    }
    if (accept)
    {
      UseCorona = true;
    }
  }
  if (strcmp(str,"CoronaSystem")==0)
  {
    in_stream >> str;
    accept = true;
    if (strcmp(str,"xrb")==0)
    {
      Corona.system = 0;
    }
    else
    {
      accept = false;
    }
    if (accept)
    {
      UseCorona = true;
    }
  }
  if (strcmp(str,"CoronaSource")==0)
  {
    in_stream >> str;
    Corona.source = (strcmp(str,"yes")==0);
    if (Corona.source || strcmp(str,"no")==0)
    {
      UseCorona = true;
      ReadCoronaSource = true;
      accept = true;
    }
  }
  if (strcmp(str,"CoronaCompton")==0)
  {
    in_stream >> str;
    Corona.compton = (strcmp(str,"yes")==0);
    if (Corona.compton || strcmp(str,"no")==0)
    {
      UseCorona = true;
      accept = true;
    }
  }
  if (strcmp(str,"CoronaPosition")==0)
  {
    for (j=0; j<=2; j++)
    {
      in_stream >> str;
      Corona.pos[j] = atof(str);
    }
    UseCorona = true;
    accept = true;
  }
  if (strcmp(str,"CoronaRadius")==0)
  {
    in_stream >> str;
    Corona.radius = atof(str);
    UseCorona = true;
    accept = true;
  }
  if (strcmp(str,"CoronaHeight")==0)
  {
    in_stream >> str;
    Corona.height = atof(str);
    UseCorona = true;
    accept = true;
  }
  if (strcmp(str,"CoronaTau")==0)
  {
    in_stream >> str;
    Corona.tau = atof(str);
    UseCorona = true;
    accept = true;
  }
  if (strcmp(str,"CoronaDensity")==0)
  {
    in_stream >> str;
    Corona.density = atof(str);
    UseCorona = true;
    accept = true;
  }
  if (strcmp(str,"CoronaTemp")==0)
  {
    in_stream >> str;
    Corona.kTe = atof(str);
    UseCorona = true;
    accept = true;
  }
  if (strcmp(str,"CoronaGamma")==0)
  {
    in_stream >> str;
    Corona.gamma = atof(str);
    UseCorona = true;
    accept = true;
  }
  if (strcmp(str,"CoronaEnergyRange")==0)
  {
    in_stream >> str;
    Corona.emin = atof(str);
    in_stream >> str;
    Corona.emax = atof(str);
    UseCorona = true;
    accept = true;
  }
  return accept;
}

bool Model::ParseDiskOption(std::istream& in_stream,char str[200],
                            bool& ReadSeedSource,bool& ReadSeedAngular,
                            bool& ReadDiskReflection,bool& ReadDiskReflectRin,
                            bool& ReadDiskReflectRout)
{
  bool accept = false;
  if (strcmp(str,"XRBSeedSource")==0)
  {
    in_stream >> str;
    SeedDisk.source = (strcmp(str,"yes")==0);
    if (SeedDisk.source || strcmp(str,"no")==0)
    {
      ReadSeedSource = true;
      accept = true;
    }
  }
  if (strcmp(str,"XRBSeedSpectrum")==0)
  {
    in_stream >> str;
    accept = true;
    if (strcmp(str,"powerlaw")==0)
    {
      SeedDisk.spectrum = 0;
    }
    else if (strcmp(str,"blackbody")==0)
    {
      SeedDisk.spectrum = 1;
    }
    else
    {
      accept = false;
    }
  }
  if (strcmp(str,"XRBSeedAngular")==0)
  {
    in_stream >> str;
    accept = true;
    if (strcmp(str,"isotropic")==0)
    {
      SeedDisk.angular = 0;
      ReadSeedAngular = true;
    }
    else if (strcmp(str,"lambert")==0)
    {
      SeedDisk.angular = 1;
      ReadSeedAngular = true;
    }
    else
    {
      accept = false;
    }
  }
  if (strcmp(str,"XRBSeedPolarization")==0)
  {
    in_stream >> str;
    accept = true;
    if (strcmp(str,"none")==0)
    {
      SeedDisk.polarization = 0;
    }
    else if (strcmp(str,"chandrasekhar")==0)
    {
      SeedDisk.polarization = 1;
    }
    else
    {
      accept = false;
    }
  }
  if (strcmp(str,"XRBSeedRin")==0)
  {
    in_stream >> str;
    SeedDisk.rin = atof(str);
    accept = true;
  }
  if (strcmp(str,"XRBSeedRout")==0)
  {
    in_stream >> str;
    SeedDisk.rout = atof(str);
    accept = true;
  }
  if (strcmp(str,"XRBSeedGamma")==0)
  {
    in_stream >> str;
    SeedDisk.gamma = atof(str);
    accept = true;
  }
  if (strcmp(str,"XRBSeedTin")==0)
  {
    in_stream >> str;
    SeedDisk.tin = atof(str);
    accept = true;
  }
  if (strcmp(str,"XRBSeedEnergyRange")==0)
  {
    in_stream >> str;
    SeedDisk.emin = atof(str);
    in_stream >> str;
    SeedDisk.emax = atof(str);
    accept = true;
  }
  if (strcmp(str,"DiskReflection")==0)
  {
    in_stream >> str;
    SeedDisk.reflection = (strcmp(str,"yes")==0);
    if (SeedDisk.reflection || strcmp(str,"no")==0)
    {
      ReadDiskReflection = true;
      accept = true;
    }
  }
  if (strcmp(str,"DiskReflectRin")==0)
  {
    in_stream >> str;
    SeedDisk.reflectRin = atof(str);
    ReadDiskReflectRin = true;
    accept = true;
  }
  if (strcmp(str,"DiskReflectRout")==0)
  {
    in_stream >> str;
    SeedDisk.reflectRout = atof(str);
    ReadDiskReflectRout = true;
    accept = true;
  }
  if (strcmp(str,"DiskReflectAlbedo")==0)
  {
    in_stream >> str;
    SeedDisk.albedo = atof(str);
    accept = true;
  }
  if (strcmp(str,"DiskReflectAngular")==0)
  {
    in_stream >> str;
    accept = true;
    if (strcmp(str,"specular")==0)
    {
      SeedDisk.reflectAngular = 0;
    }
    else if (strcmp(str,"lambert")==0)
    {
      SeedDisk.reflectAngular = 1;
    }
    else
    {
      accept = false;
    }
  }
  return accept;
}
