#ifndef disk_reflection_polarization_h
#define disk_reflection_polarization_h

#include <math.h>

struct DiskReflectionPolarization
{
  double axis[3];
  double Q,U,V;
};

inline double DiskReflectionDot(const double a[3],const double b[3])
{
  return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

// Apply the Thomson Mueller matrix in bases perpendicular to the
// incident/outgoing scattering planes. Keep the packet intensity fixed to
// preserve the existing Lambert/specular angular sampler; this corrects
// polarization transport, not the disk's angular or energy spectrum.
inline bool DiskReflectionPolarize(const double incoming[3],
                                   const double outgoing[3],
                                   const double frame1[3],
                                   const double frame2[3],
                                   double I,double Q,double U,double V,
                                   DiskReflectionPolarization& result)
{
  double norm,c,s,cos2,sin2,mu,mu2,qIn,uIn,iOut,qOut;

  result.axis[0] = incoming[1]*outgoing[2]-incoming[2]*outgoing[1];
  result.axis[1] = incoming[2]*outgoing[0]-incoming[0]*outgoing[2];
  result.axis[2] = incoming[0]*outgoing[1]-incoming[1]*outgoing[0];
  norm = sqrt(DiskReflectionDot(result.axis,result.axis));
  if (norm>1e-12)
  {
    result.axis[0] /= norm;
    result.axis[1] /= norm;
    result.axis[2] /= norm;
  }
  else
  {
    // Forward/backward scattering has no unique plane; retain the old axis.
    result.axis[0] = frame1[0];
    result.axis[1] = frame1[1];
    result.axis[2] = frame1[2];
  }

  c = DiskReflectionDot(result.axis,frame1);
  s = DiskReflectionDot(result.axis,frame2);
  cos2 = c*c-s*s;
  sin2 = 2*c*s;
  qIn = Q*cos2+U*sin2;
  uIn = -Q*sin2+U*cos2;
  mu = DiskReflectionDot(incoming,outgoing);
  if (mu>1) mu = 1;
  if (mu<-1) mu = -1;
  mu2 = mu*mu;

  iOut = 0.5*((1+mu2)*I+(1-mu2)*qIn);
  if (iOut<=1e-14*I || I<=0)
  {
    return false;
  }
  qOut = 0.5*((1-mu2)*I+(1+mu2)*qIn);
  result.Q = I*qOut/iOut;
  result.U = I*mu*uIn/iOut;
  result.V = I*mu*V/iOut;
  return true;
}

inline void DiskReflectionRotateToFrame(
    const DiskReflectionPolarization& reflected,
    const double frame1[3],const double frame2[3],
    double& Q,double& U)
{
  double c = DiskReflectionDot(reflected.axis,frame1);
  double s = DiskReflectionDot(reflected.axis,frame2);
  double norm = sqrt(c*c+s*s);
  if (norm>0)
  {
    c /= norm;
    s /= norm;
  }
  Q = reflected.Q*(c*c-s*s)-reflected.U*(2*c*s);
  U = reflected.Q*(2*c*s)+reflected.U*(c*c-s*s);
}

#endif
