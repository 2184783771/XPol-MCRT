// Implementation of the routines for a virtual detector. The object
// 'detector' adds up the Stokes vectors of the detected photons
// according to their wavelength. It then can save the stored data to
// a disk.
//
// STOKES, version 1.2, March 2014

#include <iostream>
#include <fstream>
#include <string.h>
#include <math.h>
#include "detector-stokes-v1.2.h"
#include "photon-stokes-v1.2.h"
#include "model-stokes-v1.2.h"
#include "observe-stokes-v1.2.h"

using namespace std;


// initialize detector by setting all channels to zero
Screen::Screen()
{
  I = 0;
  Q = 0;
  U = 0;
  V = 0;
  AllocRes = 0;
  AllocX = 0;
  AllocY = 0;
}


void Screen::EnsureAllocated(Model M)
{
  int i,total,res,x,y;

  res = M.res();
  x = int(M.GetResX());
  y = int(M.GetResY());
  if (res==AllocRes && x==AllocX && y==AllocY && I!=0)
  {
    return;
  }

  AllocRes = res;
  AllocX = x;
  AllocY = y;
  total = AllocRes*AllocX*AllocY;
  I = new double[total];
  Q = new double[total];
  U = new double[total];
  V = new double[total];
  for(i=0; i<total; i++)
  {
    I[i] = 0;
    Q[i] = 0;
    U[i] = 0;
    V[i] = 0;
  }
}


int Screen::Index(int channel,int x,int y)
{
  return (channel*AllocX+x)*AllocY+y;
}


double Screen::GetStokes(int n,int channel,int x,int y)
{
  int idx;

  if (I==0 || channel<0 || channel>=AllocRes ||
      x<0 || x>=AllocX || y<0 || y>=AllocY)
  {
    return 0;
  }
  idx = Index(channel,x,y);
  switch(n)
  {
    case 0:
      return I[idx];
    case 1:
      return Q[idx];
    case 2:
      return U[idx];
    default:
      return V[idx];
  }
}


// register a photon according to its wavelength and add its
// Stokes parameters
void Screen::RegisterPhoton(Photon p, Model M, double x, double y)
{
  int i,j=0,k=0,idx;
  double step,MaxExp,ResX,ResY;
  bool flag;

  EnsureAllocated(M);

  // get the maximum expansion of the XY plane and resolution degree
  MaxExp = M.GetMaxExp();
  ResX = M.GetResX();
  ResY = M.GetResY();

  // compute the wavelength space covered by one channel
  step = double(M.LambdaRange())/M.res();

  // search the correct spatial bin of the XY plane
  double X=(-MaxExp)/2;
  double Y=(-MaxExp)/2;

  while ((x>X) && (j<ResX)) 
  {
	X=X+(MaxExp/ResX);
	j=j+1;
  }
  if (j>0) j--;


  while ((y>Y) && (k<ResY))
  {
	Y=Y+(MaxExp/ResY);
        k=k+1;
  }
  if (k>0) k--;
  

  // search the correct channel and add the Stokes vector
  if ((p.GetLambda()>=M.LambdaMin()))
  {
    i = 0;
    flag = false;
    while ((i < M.res()) && !flag)
    {
      if ((p.GetLambda() < M.LambdaMin() + (i+1) * step) &&
	(x > (-MaxExp)/2) && 
        (x < MaxExp/2) &&
        (y > (-MaxExp)/2) &&   
        (y < MaxExp/2))
      {
        flag = true;
        idx = Index(i,j,k);
        I[idx] = I[idx]+p.GetI()*p.GetW();
        Q[idx] = Q[idx]+p.GetQ()*p.GetW();
        U[idx] = U[idx]+p.GetU()*p.GetW();
        V[idx] = V[idx]+p.GetV()*p.GetW();
      }
      i++;
    }
  }
}


// save the data of the detector to the disk
void Screen::SaveCurve(ofstream out_stream, int n, Model M)
{
  int i,j,k;
  double ResX,ResY;

  ResX = M.GetResX();
  ResY = M.GetResY();
  
  for (j=0; j<ResX; j++)
  {
	for (k=0; k<ResY; k++)
	{
  		for (i = 0; i < M.res(); i++)
  		{
			// the Stokes parameter to be saved is denoted by n
    			switch(n)
    			{
      				case 0:
        			  out_stream << GetStokes(0,i,j,k) << " ";
        			  break;
      				case 1:
        			  out_stream << GetStokes(1,i,j,k) << " ";
      				  break;
      				case 2:
      				  out_stream << GetStokes(2,i,j,k) << " ";
      				  break;
      				case 3:
        			  out_stream << GetStokes(3,i,j,k) << " ";
      				  break;
			}
    		}
	}
  }
}
