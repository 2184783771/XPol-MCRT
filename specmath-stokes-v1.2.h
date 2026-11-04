// This is the header file listing the definitions of all special
// mathematical routines needed. These include operations for the
// algebra of vectors, matrices, and complex numbers, as well as
// coordinate transformations, hyperbolic trigonometry, distance
// computations to various geometric  bodies, and solving algebraic
// equations of third and fourth degree. The algorithms for solving
// higher degree algebraic equations were taken from Bronstein &
// Semendjajew (1996).
//
// Reference:
// Bronstein I.N., Semendjajew K.A. 1996, Teubner Leipzig
//
// STOKES, version 1.2, March 2014


#ifndef specmath_h
#define specmath_h


// definition of pi
const double pi = 3.1415926535898;

// definition of the velocity of light in km/s
const double vc = 299792;


// type definition for complex numbers
struct complex
{
  double re;
  double im;
};


// definition of the real unit
const complex U0 = {1,0};


// limit of the imaginary part for 'real' numbers
const double ComplexLimit = 1e-4;


// return the hyperbolic arcus sine of y
double arsinh(double y);


// return the hyperbolic arcus cosine of y, if y < 1 return zero
double arcosh(double y);


// set x1, x2, and x3 to the solution of a cubic equation in x given 
// by the constants a, b and c
void SolveCubic(double a,double b,double c,complex& x1,
                complex& x2,complex& x3);


// set x1, x2, x3, and x4 to the solution of a 4th degree equation 
// in x given by the constants a, b, c, and d
void SolveFourth(double a,double b,double c,double d,complex& x1,
                 complex& x2,complex& x3,complex& x4);


// converts the spherical coordinates of a point on the unity sphere
// into Cartesian coordinates
void CartCoord(double Cartesian[],double phi,double theta);


// converts the Cartesian coordinates of a point on the unity sphere
// into spherical coordinates phi and theta
void SphericalCoords(double Cartesian[],double& phi,double& theta);


// if a trajectory through the starting point given by pos[] and in the
// direction given by direct[] hits the torus defined by R, a, and b, 
// then set t to the distance to the torus (in pc) 
void TorusCoinc(double pos[],double direct[],double& t,double a,
                double b,double R,double offset);

	
// if a trajectory through the starting point given by pos[] and in the
// direction given by direct[] hits the cylinder defined by R, a, and b, 
// then set t to the distance to the cylinder (in pc) 
void CylinderCoinc(double pos[],double direct[],double& t,
                   double r,double h,double offset);


// if a trajectory through the starting point given by pos[] and in the
// direction given by direct[] hits the ring defined by a, b, and h, 
// then set t to the distance to the ring (in pc) 
void RingCoinc(double pos[],double direct[],double& t,
               double a,double b,double h);


// if a trajectory through the starting point given by pos[] and in the
// direction given by direct[] hits the sphere double segment defined by 
// CloudPos[], r, and theta then set t to the distance to the this
// segment (in pc) 
void SphereSegCoinc(double pos[],double direct[],double CloudPos[],
                    double& t,double r,double theta,bool polar);


// if a trajectory through the starting point given by pos[] and in the
// direction given by direct[] hits the double-cone defined by a, b, 
// and theta then set t to the distance to the double-cone (in pc) 
void ConeCoinc(double pos[],double direct[],double a,double b,
               double theta,double& t);


// initialize a complex number by real part a and imaginary part b
complex C(double a,double b);


// return the real part of the complex number A
double CReal(complex A);


// return the imaginary part of the complex number A
double CImag(complex A);


// return the conjugate complex of the complex number A
complex CConj(complex A);


// return the absolute value of the complex number A
double CAbs(complex A);


// addition of two complex number A and B
complex operator +(complex A,complex B);


// addition of a real number a with a complex number B
complex operator +(double a,complex B);


// addition of a complex number A with a real number b
complex operator +(complex A,double b);


// change sign of a complex number A
complex operator -(complex A);


// subtraction of a complex number B from a complex number A
complex operator -(complex A,complex B);


// subtraction of a complex number B from a real number a
complex operator -(double a,complex B);


// subtraction of a real number b from a complex number A
complex operator -(complex A,double b);


// multiplication of two complex numbers A and B
complex operator *(complex A,complex B);


// multiplication of a real with a complex number	
complex operator *(double a,complex B);


// multiplication of a complex with a real number
complex operator *(complex A,double b);


// division of a complex number A by a complex number B
complex operator /(complex A,complex B);


// division of a real number a by a complex number B
complex operator /(double a,complex B);


// division of a complex number A by a real number b
complex operator /(complex A,double b);


// return complex square root of complex number A
complex CSqrt(complex A);


// return complex square root of real number a
complex CSqrt(double a);


// set vector a equal to vector b
void Equal(double a[],double b[]);


// set c to the vector sum of a and b
void VecSum(double c[],double a[],double b[]);


// set c to the vector difference of a and b
void VecDiff(double c[],double a[],double b[]);


// return the scalar product of vector a and vector b
double SProd(double a[],double b[]);


// set c to the s-product of the real number a and vector b
void SMult(double c[],double b,double a[]);


// set r to the product of the matrix mat with vector v
void MatVecMult(double mat[][3],double v[],double r[]);


// return the arcus tangent of y/x including the case x=y=0
double atang(double x,double y);


#endif
