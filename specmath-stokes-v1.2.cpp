// This is the implementation file of all special mathematical
// routines needed. These include operations for the algebra of
// vectors, matrices, and complex numbers, as well as coordinate
// transformations, hyperbolic trigonometry, distance computations to
// various geometric  bodies, and solving algebraic equations of third
// and fourth degree. The algorithms for solving higher degree
// algebraic equations were taken from Bronstein & Semendjajew (1996).
//
// Reference:
// Bronstein I.N., Semendjajew K.A. 1996, Teubner Leipzig
//
// STOKES, version 1.2, March 2014


#include <math.h>
#include "specmath-stokes-v1.2.h"

using namespace std;


// return the hyperbolic arcus sine of y
double arsinh(double y)
{
  return log(y+sqrt(y*y+1));
}


// return the hyperbolic arcus cosine of y, if y < 1 return zero
double arcosh(double y)
{ 
  if (y>=1)
  {
    return log(y+sqrt(y*y-1));
  }
  else
  {
    return 0;
  }
}


// set x1, x2, and x3 to the solution of a cubic equation in x given 
// by the constants a, b and c
void SolveCubic(double a,double b,double c,
		complex& x1,complex& x2,complex& x3)
{
  double p,q,D,N,beta;

  p = b/3-a*a/9;
  q = a*a*a/27-a*b/6+c/2;
  D = p*p*p+q*q;
  N = q/fabs(q)*sqrt(fabs(p));
  if (q==0)
  {
    x1 = C(0,0);
    x2 = CSqrt(-3*p);
    x3 = -CSqrt(-3*p);
  }
  else
  {
    if (p<0)
    {
      if (D<=0)
    	{
        beta = 1.0/3*acos(q/(N*N*N));
        x1 = C(-2*N*cos(beta),0);
        x2 = C(2*N*cos(beta+pi/3),0);
        x3 = C(2*N*cos(beta-pi/3),0);
      }
      else
      {
	beta = 1.0/3*arcosh(q/(N*N*N));
     	x1 = C(-2*N*cosh(beta),0);
     	x2 = C(N*cosh(beta),N*sqrt(3.0)*sinh(beta));
     	x3 = C(N*cosh(beta),-N*sqrt(3.0)*sinh(beta));
      }
    }
    else
    {
      beta = 1.0/3*arsinh(q/(N*N*N));
      x1 = C(-2*N*sinh(beta),0);
      x2 = C(N*sinh(beta),N*sqrt(3.0)*cosh(beta));
      x3 = C(N*sinh(beta),-N*sqrt(3.0)*cosh(beta));
    }
  }
  x1 = x1-a/3;
  x2 = x2-a/3;
  x3 = x3-a/3;
}


// set x1, x2, x3, and x4 to the solution of a 4th degree equation 
// in x given by the constants a, b, c, and d
void SolveFourth(double a,double b,double c,double d,
		complex& x1,complex& x2,complex& x3,complex& x4)
{
  double p,q,r;
  complex alpha,beta,gamma,u,v,w;

  p = b-3.0/8*a*a;
  q = 1.0/8*a*a*a-1.0/2*a*b+c;
  r = d-1.0/4*a*c+1.0/16*a*a*b-3.0/256*a*a*a*a;
  SolveCubic(2*p,p*p-4*r,-q*q,alpha,beta,gamma);
  u = CSqrt(alpha);
  v = CSqrt(beta);
  w = CSqrt(gamma);
  if (CReal(u*v*w)*q<0)
  {
    u=-u;
  }
  if (CReal(u*v*w)*q<0)
  {
    v=-v;
  }
  if (CReal(u*v*w)*q<0)
  {
    w=-w;
  }
  x1 = 1.0/2*(u+v-w)-a/4;
  x2 = 1.0/2*(u-v+w)-a/4;
  x3 = 1.0/2*(-u+v+w)-a/4;
  x4 = 1.0/2*(-u-v-w)-a/4;
}


// converts the spherical coordinates of a point on the unity sphere
// into Cartesian coordinates
void CartCoord(double Cartesian[],double phi,double theta)
{
  Cartesian[0] = cos(phi)*sin(theta);
  Cartesian[1] = sin(phi)*sin(theta);
  Cartesian[2] = cos(theta);
}


// converts the Cartesian coordinates of a point on the unity sphere
// into spherical coordinates phi and theta
void SphericalCoords(double Cartesian[],double& phi,double& theta)
{
  double cs;

  if (fabs(Cartesian[2])>=1)
  {
    if (Cartesian[2]>0)
    {
      theta = 0;
    }
    else
    {
      theta = pi;
    }
    phi = 0;
  }
  else
  {
    theta = acos(Cartesian[2]);
    cs = Cartesian[0]/sin(theta);
    if (fabs(cs)<=1)
    {
      if (Cartesian[1]>0)
      {
        phi=acos(cs);
      }
      else
      {
        phi=2*pi-acos(cs);
      }
    }
    else
    {
      if (cs>0)
      {
        phi=0;
      }
      else
      {
        phi=pi;
      }
    }
  }
}


// if a trajectory through the starting point given by pos[] and in the
// direction given by direct[] hits the torus defined by R, a, and b, 
// then set t to the distance to the torus (in pc) 
void TorusCoinc(double pos[],double direct[],double& t,double a,
                double b,double R,double offset)
{
  int i;
  double a0,a1,a2,a3,c1,c2,c3,c4;
  complex x[4];

  c1 = a*a*(direct[0]*direct[0]+direct[1]*direct[1])+b*b*direct[2]*direct[2];
  c2 = a*a*(direct[0]*pos[0]+direct[1]*pos[1])+b*b*direct[2]*(pos[2]-offset);
  c3 = a*a*(R*R+pos[0]*pos[0]+pos[1]*pos[1])+b*b*(pos[2]-offset)*(pos[2]-offset)-a*a*b*b;
  c4 = 4*a*a*a*a*R*R;
  a3 = 4*c2/c1;
  a2 = (2*c1*c3+4*c2*c2-c4*(direct[0]*direct[0]+direct[1]*direct[1]))/(c1*c1);
  a1 = (4*c2*c3-2*c4*(direct[0]*pos[0]+direct[1]*pos[1]))/(c1*c1);
  a0 = (c3*c3-c4*(pos[0]*pos[0]+pos[1]*pos[1]))/(c1*c1);
  SolveFourth(a3,a2,a1,a0,x[0],x[1],x[2],x[3]);
  for (i=0; i<4; i++)
  {
    if ((fabs(CImag(x[i]))<ComplexLimit) && (CReal(x[i])>0) && (CReal(x[i])<t))
    {
      t = CReal(x[i]);
    }
  }
}


// if a trajectory through the starting point given by pos[] and in the
// direction given by direct[] hits the cylinder defined by R, a, and b, 
// then set t to the distance to the cylinder (in pc) 
void CylinderCoinc(double pos[],double direct[],double& t,
		   double r,double h,double offset)
{
  double p,q,z,t0,t1,t2;

  t0 = t;
  if ((direct[0]!=0)||(direct[1]!=0))
  {
    p = (pos[0]*direct[0]+pos[1]*direct[1])/
        (direct[0]*direct[0]+direct[1]*direct[1]);
    q = (pos[0]*pos[0]+pos[1]*pos[1]-r*r)/
        (direct[0]*direct[0]+direct[1]*direct[1]);
    t2 = p*p-q;
    if (t2>=0)
    {
      t2 = sqrt(t2);
      t1 = -p-t2;
      z = pos[2]-offset+t1*direct[2];
      if ((t1>=0)&&(fabs(z)<=h))
      {
        t0 = t1;
      }
      else
      {
        t2 = t1 + 2*t2;
        z = pos[2]-offset+t2*direct[2];
        if ((t2>=0)&&(fabs(z)<=h))
        {
          t0 = t2;
        }
      }
    }
  }
  else
  {
    if ((sqrt(pos[0]*pos[0]+pos[1]*pos[1])==r)&&(fabs(pos[2]-offset)<=h))
    {
      t0 = 0;
    }
  }
  if (t0<t)
  {
    t = t0;
  }
}


// if a trajectory through the starting point given by pos[] and in the
// direction given by direct[] hits the ring defined by a, b, and h, 
// then set t to the distance to the ring (in pc) 
void RingCoinc(double pos[],double direct[],double& t,
               double a,double b,double h)
{
  double r,t0;

  if ((direct[2]!=0))
  {
    t0 = (h-pos[2])/direct[2];
    if (t0>=0)
    {
      r = sqrt(pow(pos[0]+t0*direct[0],2) +
           pow(pos[1]+t0*direct[1],2));
      if ((r>=a)&&(r<=b))
      {
        if (t0<t)
        {
          t = t0;
        }
      }
    }
  }
}


// if a trajectory through the starting point given by pos[] and in
// the direction given by direct[] hits the sphere double segment
// defined by  CloudPos[], r, and theta then set t to the distance to
// the this segment (in pc). Theta may define either the opening-angle
// of a double-cone or of a flared disk.
void SphereSegCoinc(double pos[],double direct[],double CloudPos[],
		    double& t,double r,double theta,bool polar)
{
  double t0,t1,t2,tn,con[3];

  VecDiff(con,pos,CloudPos);
  t1 = SProd(con,direct);
  t2 = SProd(con,con);
  t2 = t1*t1 - t2 + r*r;
  if (t2>=0)
  {
    t0 = t;
    t1 = -t1 - sqrt(t2);
    if (t1>=0)
    {
      t0 = t1;
    }
    else
    {
      t2 = t1 + 2*sqrt(t2);
      if (t2>=0)
      {
        t0 = t2;
      }
    }
    if (t0<t)
    {
      if (theta<pi/2)
      {
        SMult(con,t0,direct);
        VecSum(con,pos,con);
	if (polar)
	{
          tn = tan(pi/2-theta);
          if(tn*tn*(con[0]*con[0]+con[1]*con[1]) < con[2]*con[2])
          {
            t = t0;
          }
	}
	else
	{
          tn = tan(theta);
          if(tn*tn*(con[0]*con[0]+con[1]*con[1]) > con[2]*con[2])
          {
            t = t0;
          }
	}
      }
      else
      {
        t = t0;
      }
    }
  }
}


// if a trajectory through the starting point given by pos[] and in the
// direction given by direct[] hits the double-cone defined by a, b, 
// and theta then set t to the distance to the this double-cone (in pc) 
void ConeCoinc(double pos[],double direct[],double a,double b,
	       double theta,double& t)
{
  double d,p,q,l,tn,t0,t1,con[3];

  t0 = t;
  tn = tan(pi/2-theta);
  d = direct[2]*direct[2]-tn*tn*(direct[0]*direct[0]+direct[1]*direct[1]);
  if (d!=0)
  {
    p = (pos[2]*direct[2]-tn*tn*(pos[0]*direct[0]+pos[1]*direct[1]))/d;
    q = (pos[2]*pos[2]-tn*tn*(pos[0]*pos[0]+pos[1]*pos[1])) / d;
    q = p*p-q;
    if (q>=0)
    {
      q = sqrt(q);
      t1 = -p - q;
      SMult(con,t1,direct);
      VecSum(con,pos,con);
      l = SProd(con,con);
      if ((t1>=0)&&(l>=a*a)&&(l<=b*b))
      {
        t0 = t1;
      }
      else
      {
        t1 = t1 + 2*q;
        SMult(con,t1,direct);
        VecSum(con,pos,con);
        l = SProd(con,con);
        if ((t1>=0)&&(l>=a*a)&&(l<=b*b))
        {
          t0 = t1;
        }
      }
    }
  }
  else
  {
    l = SProd(pos,pos);
    if ((tn*tn*(pos[0]*pos[0]+pos[1]*pos[1])==pos[2]*pos[2])&&
        (l>=a*a)&&(l<=b*b))
    {
      t0 = 0;
    }
  }
  if (t0<t)
  {
    t = t0;
  }
}


// initialize a complex number
complex C(double a,double b)
{
  complex res;

  res.re=a;
  res.im=b;
  return res;
}


// return the real part of the complex number A
double CReal(complex A)
{
  return A.re;
}


// return the imaginary part of the complex number A
double CImag(complex A)
{
  return A.im;
}


// return the conjugate complex of the complex number A
complex CConj(complex A)
{
  complex res;

  res.re = A.re;
  res.im = -A.im;
  return res;
}


// return the absolute value of the complex number A
double CAbs(complex A)
{
  return sqrt((A*CConj(A)).re);
}


// addition of two complex number A and B
complex operator +(complex A,complex B)
{
  complex res;
  
  res.re = A.re + B.re;
  res.im = A.im + B.im;
  return res;
}


// addition of a real number a with a complex number B
complex operator +(double a,complex B)
{
  return a*U0+B;
}


// addition of a complex number A with a real number b
complex operator +(complex A,double b)
{
  return A+b*U0;  
}


// change sign of a complex number A
complex operator -(complex A)
{
  return -1*A;
}


// subtraction of a complex number B from a complex number A
complex operator -(complex A,complex B)
{
  return A+(-B);
}


// subtraction of a complex number B from a real number a
complex operator -(double a,complex B)
{
  return a+(-B);
}


// subtraction of a real number b from a complex number A
complex operator -(complex A,double b)
{
  return A+(-b);
}


// multiplication of two complex numbers A and B
complex operator *(complex A,complex B)
{
  complex res;

  res.re = A.re*B.re - A.im*B.im;
  res.im = A.re*B.im + B.re*A.im;
  return res;
}


// multiplication of a real with a complex number    
complex operator *(double a,complex B)
{
  complex res;
  
  res.re = a*B.re;
  res.im = a*B.im;
  return res;
}


// multiplication of a complex with a real number
complex operator *(complex A,double b)
{
  return b*A;
}


// division of a complex number A by a complex number B
complex operator /(complex A,complex B)
{
  return 1/pow(CAbs(B),2)*A*CConj(B);
}


// division of a real number a by a complex number B
complex operator /(double a,complex B)
{
  return a*U0/B;
}


// division of a complex number A by a real number b
complex operator /(complex A,double b)
{
  return 1/b*A;
}


// return complex square root of complex number A
complex CSqrt(complex A)
{
  double r,phi;

  r = sqrt(CAbs(A));
  phi = 0.5*atang(A.im,A.re);
  return C(r*cos(phi),r*sin(phi));
}


// return complex square root of real number a
complex CSqrt(double a)
{
  return CSqrt(C(a,0));
}


// set vector a equal to vector b
void Equal(double a[],double b[])
{
  int i;

  for(i=0; i<=2; i++)
    a[i] = b[i];
}


// set c to the vector sum of a and b
void VecSum(double c[],double a[],double b[])
{
  int i;

  for(i=0; i<=2; i++)
    c[i] = a[i] + b[i];
}


// set c to the vector difference of a and b
void VecDiff(double c[],double a[],double b[])
{
  int i;

  for(i=0; i<=2; i++)
    c[i] = a[i] - b[i];
}


// return the scalar product of vector a and vector b
double SProd(double a[],double b[])
{
  int i;
  double sp;

  sp = 0;
  for(i=0; i<=2; i++)
    sp = sp + a[i]*b[i];

  return sp;
}


// set c to the s-product of the real number a and vector b
void SMult(double c[],double b,double a[])
{
  int i;

  for(i=0; i<=2; i++)
    c[i] = b*a[i];
}


// set r to the product of the matrix mat with vector v
void MatVecMult(double mat[][3],double v[],double r[])
{
  int i,j;

  for(i = 0; i <= 2; i++) r[i] = 0;
  for(i = 0; i <= 2; i++)
  {
    for(j = 0; j <= 2; j++) r[i] = r[i]+mat[i][j]*v[j];
  }
}


// return the arcus tangent of y/x including the case x=y=0
double atang(double x,double y)
{
  if ((x==0)&&(y==0))
  {
    return 0;
  }
  else
  {
    return atan2(x,y);
  }
}
