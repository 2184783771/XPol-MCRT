// Implementation of the routines for the virtual detector system. The
// object 'observe' manages a spherical web of detectors. It contains
// routines to register photons escaping from the model region and to
// save the acquired polarimetric data of each detector to a disk.
//
// STOKES, version 1.2, March 2014

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <string.h>
#include <math.h>
#include <time.h>
#ifdef _MSC_VER
  #include <direct.h>
#else
  #include <sys/stat.h>
#endif
#include "observe-stokes-v1.2.h"
#include "detector-stokes-v1.2.h"

using namespace std;

const double keVAngstrom = 12.3984193;


static void EnsureOutputDirectory(const char outputfile[])
{
  int i,last,length;
  char directory[300];

  last = -1;
  for(i=0; outputfile[i]!='\0'; i++)
  {
    if (outputfile[i]=='/' || outputfile[i]=='\\')
    {
      last = i;
    }
  }
  length = int(strlen(outputfile));
  if (last<0 || last>=299 || length>=300)
  {
    return;
  }

  strcpy(directory,outputfile);
  for(i=1; i<=last; i++)
  {
    if (directory[i]=='/' || directory[i]=='\\')
    {
      char separator = directory[i];
      directory[i] = '\0';
      if (!(i==2 && directory[1]==':'))
      {
#ifdef _MSC_VER
        _mkdir(directory);
#else
        mkdir(directory,0755);
#endif
      }
      directory[i] = separator;
    }
  }
  directory[last] = '\0';

#ifdef _MSC_VER
  _mkdir(directory);
#else
  mkdir(directory,0755);
#endif
}



// initialize detectors by setting all channels to zero
Monitor::Monitor()
{
  Col = new Column[ScreenNumThetaMax];
  DirectCol = new Column[ScreenNumThetaMax];
  ReflectedCol = new Column[ScreenNumThetaMax];
  CoronaScatteredCol = new Column[ScreenNumThetaMax];
}


// register photon by the detector in its flight direction
// and save its wavelength and Stokes vector
void Monitor::RegisterPhoton(Photon p,Model M)
{
  RegisterPhotonToColumn(p,M,Col);
  if (p.IsReflectedComponent())
  {
    RegisterPhotonToColumn(p,M,ReflectedCol);
  }
  else
  {
    RegisterPhotonToColumn(p,M,DirectCol);
  }
  if (p.GetCoronaScatterNum()>0)
  {
    RegisterPhotonToColumn(p,M,CoronaScatteredCol);
  }
}


void Monitor::RegisterPhotonToColumn(Photon p,Model M,ColumnPtr Data)
{
  int i,j,ThetaNum,PhiNum;
  double Coords[3],phi,theta,Position[3],Pos[3],phi2,theta2,phi0,theta0,ro,MaxExp,ResX,ResY,x,y;
  const double pi=3.14159265358979323846264338327950;

  MaxExp = M.GetMaxExp();
  ResX = M.GetResX();
  ResY = M.GetResY();

  // get the number of detectors in the theta and in the phi direction
  PhiNum = M.ScreenNumPhi();
  ThetaNum = M.ScreenNumTheta();

  // compute angular difference between neighboring viewing angles
  PhiStep = 2*pi / PhiNum;
  ThetaStep = pi / ThetaNum;

  // get the photon's flight direction
  Coords[0] = p.GetDirect1();
  Coords[1] = p.GetDirect2();
  Coords[2] = p.GetDirect3();
  SphericalCoords(Coords,phi,theta);

  // get the photon's flight position
  Position[0] = p.GetPos1();
  Position[1] = p.GetPos2();
  Position[2] = p.GetPos3();

  // calculation of ro linked with the plane projection for imaging
  ro = sqrt((Position[0]*Position[0])+(Position[1]*Position[1])+(Position[2]*Position[2]));

  // safety condition for ro=0
  if (ro!=0)
  {
	Pos[0] = Position[0]/ro;
  	Pos[1] = Position[1]/ro;
  	Pos[2] = Position[2]/ro;
  }
  else
  {
	Pos[0] = Position[0]/0.000001;
  	Pos[1] = Position[1]/0.000001;
  	Pos[2] = Position[2]/0.000001;
  }  
  SphericalCoords(Pos,phi0,theta0);


  // if the model space is symmetric with respect to the xy-plane
  // combine the results of the two half-spaces
  if (M.GetPlaneSymmetry())
  {
    if (Coords[2] < 0)
    {
	p.TurnDirection(0,pi);
	Coords[2] = -Coords[2];
    }

    // find the correct detector in the theta direction
    i = int(floor((1-Coords[2])*ThetaNum));

  }

  // else treat all theta directions between 0 and pi
  else
  {

    // find the correct detector in the theta direction of the photon
    i = int(floor(0.5*(1-Coords[2])*ThetaNum));

  }

  // correct i if it is off-limits due to numerical inaccuracies
  if (i<0)
  {
     i=0;
  }
  if (i>=ThetaNum)
  {
     i=ThetaNum-1;
  }

  // adjust the photon's polarization plane to the detector
  p.AdjustFrame(phi,theta);

  // Rotate the model detector basis into the celestial sky basis.  The
  // position angle is measured from celestial north toward east.  A value
  // of zero preserves the historical STOKES output exactly.
  p.TurnPolarPlane(pi/180.0*M.GetSkyAxisPA());

  // find the correct detector in the phi direction
  j = int(floor(phi/PhiStep));
  if (j<0)
  {
    j=0;
  }
  if (j>=PhiNum)
  {
     j=PhiNum-1;
  }

   // calculation of the detector's angles
  if (M.GetPlaneSymmetry())
  {
    theta2=acos(1-(i+0.5)/ThetaNum);
  }
  else
  {
    theta2=acos(1-2*(i+0.5)/ThetaNum);
  }
  phi2=2*pi*((j+0.5)/PhiNum);


  // calculation of XY plane projection's coordinates 
  x=-ro*sin(theta0)*sin(phi2-phi0);
  y=ro*(cos(theta2)*sin(theta0)*cos(phi2-phi0)-cos(theta0)*sin(theta2));

  // register the photon
  if (M.GetAxisSymmetry())
  {
    Data[i].Row[0].RegisterPhoton(p,M,x,y);
  }
  else
  {
    Data[i].Row[j].RegisterPhoton(p,M,x,y);
  }

}


// save the data of all detectors to a disk
void Monitor::SaveResults(long int Photons,Model M)
{
  if (M.HasInclination())
  {
    SaveBandSummary(Photons,M,Col,"");
    if (M.ShouldSaveComponents())
    {
      SaveBandSummary(Photons,M,DirectCol,"_direct");
      SaveBandSummary(Photons,M,ReflectedCol,"_reflected");
      SaveBandSummary(Photons,M,CoronaScatteredCol,"_corona_scattered");
    }
  }
  if (M.ShouldSaveStokesFiles())
  {
    SaveDetectorSet(Photons,M,Col,"");
    if (M.ShouldSaveComponents())
    {
      SaveDetectorSet(Photons,M,DirectCol,"_direct");
      SaveDetectorSet(Photons,M,ReflectedCol,"_reflected");
      SaveDetectorSet(Photons,M,CoronaScatteredCol,"_corona_scattered");
    }
  }
  if (M.ShouldSaveMap())
  {
    SaveMapSet(Photons,M,Col,"_total");
    if (M.ShouldSaveComponents())
    {
      SaveMapSet(Photons,M,DirectCol,"_direct");
      SaveMapSet(Photons,M,ReflectedCol,"_reflected");
      SaveMapSet(Photons,M,CoronaScatteredCol,"_corona_scattered");
    }
  }
}


void Monitor::SaveBandSummary(long int Photons,Model M,ColumnPtr Data,
                              const char suffix[])
{
  const double pi=3.14159265358979323846264338327950;
  int phiNum = M.GetAxisSymmetry() ? 1 : M.ScreenNumPhi();
  int thetaIndex = 0;
  int phiIndex = 0;
  double best = 1.0e99;
  double targetTheta = M.GetInclinationDeg()*pi/180.0;
  double energy,lambda,theta,matchedTheta=0.0,mu;
  const double bandMin[4]={2.0,3.0,5.0,2.0};
  const double bandMax[4]={3.0,5.0,8.0,8.0};
  double I[4]={0,0,0,0},Q[4]={0,0,0,0};
  double U[4]={0,0,0,0},V[4]={0,0,0,0};
  int used[4]={0,0,0,0};
  char path[300];
  ofstream out;
  for (int t=0; t<M.ScreenNumTheta(); ++t)
  {
    if (M.GetPlaneSymmetry())
      mu = 1.0-(t+0.5)/M.ScreenNumTheta();
    else
      mu = 1.0-2.0*(t+0.5)/M.ScreenNumTheta();
    if (mu>1.0) mu=1.0;
    if (mu<-1.0) mu=-1.0;
    theta=acos(mu);
    if (fabs(theta-targetTheta)<best)
    {
      best=fabs(theta-targetTheta);
      thetaIndex=t;
      matchedTheta=theta;
    }
  }
  for (int channel=0; channel<M.res(); ++channel)
  {
    lambda=M.LambdaMin()+(channel+0.5)*M.LambdaRange()/M.res();
    energy=lambda>0.0 ? keVAngstrom/lambda : 0.0;
    for (int band=0; band<4; ++band)
    {
      bool inBand=(energy>=bandMin[band] &&
                   (energy<bandMax[band] || band==3 &&
                    energy<=bandMax[band]));
      if (inBand)
      {
        for (int phi=0; phi<phiNum; ++phi)
        {
          I[band] += Data[thetaIndex].Row[phi].GetStokes(0,channel,0,0);
          Q[band] += Data[thetaIndex].Row[phi].GetStokes(1,channel,0,0);
          U[band] += Data[thetaIndex].Row[phi].GetStokes(2,channel,0,0);
          V[band] += Data[thetaIndex].Row[phi].GetStokes(3,channel,0,0);
        }
        ++used[band];
      }
    }
  }
  M.GetOutputFile(path);
  strcat(path,suffix);
  strcat(path,"_energy_band_components.csv");
  EnsureOutputDirectory(path);
  out.open(path,ios::out);
  out << "component,energy_min_keV,energy_max_keV,inclination_input_deg,inclination_view_deg,photons_simulated,I,Q,U,V,PD_percent,PD_error_percent,PA_deg,PA_error_deg,spectral_bins_used\n";
  out << setprecision(10);
  for (int band=0; band<4; ++band)
  {
    double pd=0.0,pa=0.0,pdErr=0.0,paErr=0.0;
    if (I[band]>0.0 && used[band]>0)
    {
      double polarized=sqrt(Q[band]*Q[band]+U[band]*U[band]);
      pd=polarized/I[band];
      pa=0.5*atan2(U[band],Q[band])*180.0/pi;
      if (polarized>0.0)
      {
        pdErr=sqrt(polarized*polarized/(I[band]*I[band]*I[band]));
        paErr=0.5*sqrt(I[band])/polarized*180.0/pi;
      }
    }
    out << (suffix[0]=='\0' ? "total" : suffix+1) << ","
        << bandMin[band] << "," << bandMax[band] << ","
        << M.GetInclinationDeg() << "," << matchedTheta*180.0/pi << ","
        << Photons << "," << I[band] << "," << Q[band] << ","
        << U[band] << "," << V[band] << "," << 100.0*pd << ","
        << 100.0*pdErr << "," << pa << "," << paErr << ","
        << used[band] << "\n";
  }
  out.close();
}


void Monitor::SaveMapSet(long int Photons,Model M,ColumnPtr Data,
                         const char suffix[])
{
  int thetaIndex = M.GetMapThetaIndex();
  int PhiNum = M.GetAxisSymmetry() ? 1 : M.ScreenNumPhi();
  int phiIndex = 0;
  int ResX = int(M.GetResX());
  int ResY = int(M.GetResY());
  int UsedChannels = 0;
  double EnergyMin = M.GetMapEnergyMin();
  double EnergyMax = M.GetMapEnergyMax();
  double MaxExp = M.GetMaxExp();
  double I,Q,U,V,P,PA,Energy,Lambda,ThetaDeg,Mu;
  double HalfExp = MaxExp/2.0;
  char outputfile[300];
  ofstream out_stream;
  vector<int> Channels;

  if (thetaIndex < 0) thetaIndex = 0;
  if (thetaIndex >= M.ScreenNumTheta()) thetaIndex = M.ScreenNumTheta()-1;

  if (M.GetPlaneSymmetry())
  {
    Mu = 1.0 - (thetaIndex+0.5)/M.ScreenNumTheta();
  }
  else
  {
    Mu = 1.0 - 2.0*(thetaIndex+0.5)/M.ScreenNumTheta();
  }
  if (Mu>1.0) Mu = 1.0;
  if (Mu<-1.0) Mu = -1.0;
  ThetaDeg = 180.0/pi*acos(Mu);

  for (int channel=0; channel<M.res(); channel++)
  {
    Lambda = M.LambdaMin() + (channel+0.5)*M.LambdaRange()/M.res();
    Energy = (Lambda>0.0) ? keVAngstrom/Lambda : 0.0;
    if (Energy>=EnergyMin && Energy<EnergyMax)
    {
      Channels.push_back(channel);
    }
  }

  M.GetOutputFile(outputfile);
  strcat(outputfile,suffix);
  strcat(outputfile,"_map.csv");
  EnsureOutputDirectory(outputfile);
  out_stream.open(outputfile,ios::out);
  ostringstream rows;
  rows << setprecision(8);

  rows << "# XPol-MCRT band-integrated spatial map\n";
  rows << "# photons_saved," << Photons << "\n";
  rows << "# photons_total," << M.PhotonNum() << "\n";
  rows << "# theta_index," << thetaIndex << "\n";
  rows << "# theta_deg," << ThetaDeg << "\n";
  rows << "# phi_index," << phiIndex << "\n";
  rows << "# energy_min_keV," << EnergyMin << "\n";
  rows << "# energy_max_keV," << EnergyMax << "\n";
  rows << "# res_x," << ResX << "\n";
  rows << "# res_y," << ResY << "\n";
  rows << "# max_exp," << MaxExp << "\n";
  rows << "x_index,y_index,x,y,I,Q,U,V,PD,PA_deg,used_channels\n";

  for (int x=0; x<ResX; x++)
  {
    for (int y=0; y<ResY; y++)
    {
      I = Q = U = V = 0.0;
      UsedChannels = 0;
      for (size_t index=0; index<Channels.size(); index++)
      {
        int channel = Channels[index];
        I += Data[thetaIndex].Row[phiIndex].GetStokes(0,channel,x,y);
        Q += Data[thetaIndex].Row[phiIndex].GetStokes(1,channel,x,y);
        U += Data[thetaIndex].Row[phiIndex].GetStokes(2,channel,x,y);
        V += Data[thetaIndex].Row[phiIndex].GetStokes(3,channel,x,y);
        UsedChannels++;
      }
      if (I>0.0)
      {
        P = sqrt(Q*Q+U*U)/I;
        PA = 0.5*180.0/pi*atan2(U,Q);
      }
      else
      {
        P = 0.0;
        PA = 0.0;
      }
      rows << x << "," << y << ",";
      rows << -HalfExp + (x+0.5)*MaxExp/ResX << ",";
      rows << -HalfExp + (y+0.5)*MaxExp/ResY << ",";
      rows << I << "," << Q << "," << U << "," << V << ",";
      rows << P << "," << PA << "," << UsedChannels << "\n";
    }
  }
  out_stream << rows.str();
  out_stream.close();
}


void Monitor::SaveDetectorSet(long int Photons,Model M,ColumnPtr Data,
                              const char suffix[])
{
  int i,j,n, PhiNum;
  char outputfile[300];
  CTRegion CTR;
  BLRegion BLR;
  NLRegion NLR;
  scatterer Cloud;
  ofstream out_stream[4];
  time_t rawtime;
  double ResX = M.GetResX();
  double ResY = M.GetResY();
  const int StokesNum = M.ShouldSaveV() ? 4 : 3;
  const int OutputBufferSize = 1024*1024;
  static char OutputBuffer[4][1024*1024];

  // prepare file names for the output files, there is one file for
  // each Stokes parameter
  for(n=0; n<StokesNum; n++)
  {
    M.GetOutputFile(outputfile);
    strcat(outputfile,suffix);
    out_stream[n].rdbuf()->pubsetbuf(OutputBuffer[n], OutputBufferSize);
    switch(n)
    {
    case 0:
      strcat(outputfile,"I.dat");
      EnsureOutputDirectory(outputfile);
      out_stream[n].open(outputfile,ios::out);
      break;
    case 1:
      strcat(outputfile,"Q.dat");
      EnsureOutputDirectory(outputfile);
      out_stream[n].open(outputfile,ios::out);
      break;
    case 2:
      strcat(outputfile,"U.dat");
      EnsureOutputDirectory(outputfile);
      out_stream[n].open(outputfile,ios::out);
      break;
    default:
      strcat(outputfile,"V.dat");
      EnsureOutputDirectory(outputfile);
      out_stream[n].open(outputfile,ios::out);
      break;
    }
    out_stream[n] << setprecision(8);
    outputfile[strlen(outputfile)-4]='\0';

    // save all parameters of the model run -

    // ...model run, program version, date, time and scattering medium
    time(&rawtime);
    out_stream[n] << outputfile << " XPol-MCRT version 1.0 "
                  << "electron " << ctime(&rawtime) << "\n";

    // ... emitted spectrum:
    out_stream[n] << Photons << " " << M.PhotonNum() << " ";
    out_stream[n] << M.LambdaMin() << " ";
    out_stream[n] << M.LambdaMax() << " ";
    out_stream[n] << M.LambdaRange() << " ";

    // ...continuum emission region
    CTR = M.GetCTR();
    out_stream[n] << CTR.R << " ";
    out_stream[n] << CTR.a << " ";
    out_stream[n] << CTR.b << " ";
    out_stream[n] << CTR.c << " ";
    out_stream[n] << CTR.alpha << " ";
    out_stream[n] << CTR.vtype << " ";
    out_stream[n] << CTR.vr << " ";
    out_stream[n] << CTR.vphi << " ";
    out_stream[n] << CTR.vtheta << " ";

    // ...broad line emission region
    BLR = M.GetBLR();
    out_stream[n] << BLR.R << " ";
    out_stream[n] << BLR.a << " ";
    out_stream[n] << BLR.b << " ";
    out_stream[n] << BLR.lambda << " ";
    out_stream[n] << BLR.width << " ";
    out_stream[n] << BLR.rli << " ";
    out_stream[n] << BLR.vtype << " ";
    out_stream[n] << BLR.vr << " ";
    out_stream[n] << BLR.vphi << " ";
    out_stream[n] << BLR.vtheta << " ";


    // ...narrow line emission region
    NLR = M.GetNLR();
    out_stream[n] << NLR.a << " ";
    out_stream[n] << NLR.b << " ";
    out_stream[n] << NLR.theta << " ";
    out_stream[n] << NLR.lambda << " ";
    out_stream[n] << NLR.width << " ";
    out_stream[n] << NLR.rli << " ";
    out_stream[n] << NLR.vtype << " ";
    out_stream[n] << NLR.vr << " ";
    out_stream[n] << NLR.vphi << " ";
    out_stream[n] << NLR.vtheta << " ";


    // ...web detector
    out_stream[n] << M.ScreenNumPhi() << " ";
    out_stream[n] << M.ScreenNumTheta() << " ";
    out_stream[n] << M.res() << " ";
    out_stream[n] << M.GetPlaneSymmetry() << " ";
    out_stream[n] << M.GetAxisSymmetry() << " ";

    // ...imaging resolution
    out_stream[n] << M.GetResX() << " ";
    out_stream[n] << M.GetResY() << " ";
    out_stream[n] << M.GetMaxExp() << " ";

    // ...computational parameters
    out_stream[n] << M.GetGap() << " ";

    // save data of the scattering regions defined
    out_stream[n] << M.GetScattNum() << " ";
    for(i=0; i<M.GetScattNum(); i++)
    {
    	Cloud = M.GetCloud(i);
    	out_stream[n] << Cloud.shape << " ";
     	switch(Cloud.shape)
      	{
        	case 0:
         	 out_stream[n] << Cloud.R << " ";
         	 out_stream[n] << Cloud.a << " ";
         	 out_stream[n] << Cloud.b << " ";
         	 break;
       		case 1:
        	 out_stream[n] << Cloud.R << " ";
          	 out_stream[n] << Cloud.a << " ";
         	 out_stream[n] << Cloud.b << " ";
         	 out_stream[n] << Cloud.c << " ";
         	 break;
       		case 2:
         	 out_stream[n] << Cloud.a << " ";
         	 out_stream[n] << Cloud.b << " ";
         	 out_stream[n] << 180/pi*Cloud.theta << " ";
         	 break;
       		case 3:
          	 for(j=0; j<=2; j++) out_stream[n] << Cloud.pos[j] << " ";
          	 out_stream[n] << Cloud.a << " ";
         	 break;
       		default:
          	 out_stream[n] << Cloud.a << " ";
          	 out_stream[n] << Cloud.b << " ";
          	 out_stream[n] << 180/pi*Cloud.theta << " ";
         	 break;
      	}
      	out_stream[n] << Cloud.Density << " ";
        out_stream[n] << Cloud.vtype << " ";
        out_stream[n] << Cloud.vr << " ";
        out_stream[n] << Cloud.vphi << " ";
        out_stream[n] << Cloud.vtheta << " ";
    }

    // save tables of the Stokes parameters -
    out_stream[n] << "\n\n*****\n\n";


   // vertical: viewing direction in theta and phi
   if (M.GetAxisSymmetry())
   {
     PhiNum = 1;
   }
   else
   {
     PhiNum = M.ScreenNumPhi();
   }
   for (i=0; i<M.ScreenNumTheta(); i++)
   {
   	for (j=0; j<PhiNum; j++)
      	{
    		for (int k=0; k<ResX; k++)
    		{
			for (int l=0; l<ResY; l++)
      			{
      				// horizontal: spectroscopic channels.  Buffer each row before
              // writing to disk; this is much faster for imaging output.
              ostringstream row_stream;
              row_stream << setprecision(8);
      				for(int itmp=0; itmp < M.res(); itmp++)
        			{ 
          				switch(n)
          				{
            					case 0:
              					  row_stream << Data[i].Row[j].GetStokes(0,itmp,k,l) << " ";
            					  break;
            					case 1:
            					  row_stream << Data[i].Row[j].GetStokes(1,itmp,k,l) << " ";
            					  break;
            					case 2:
            					  row_stream << Data[i].Row[j].GetStokes(2,itmp,k,l) << " ";
            					  break;
            					default:
            					  row_stream << Data[i].Row[j].GetStokes(3,itmp,k,l) << " ";
             					  break;
         				}
        			}
              row_stream << "\n";
        			out_stream[n] << row_stream.str();
      			}
      			out_stream[n] << "\n";
    		}
	}
    }
    out_stream[n].close();
  }
}
