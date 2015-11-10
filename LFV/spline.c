#include "mynr.h"
#include "mynrutil.h"

/* Construct a cubic spline.

   Input:
   x --- vector of x values.
   y --- vector of y values.
   n --- number of data points (vector length).
   yp1 --- dy/dx(x[0]).
   ypn --- dy/dx(x[n-1]).
   Output:
   y2 --- vector of 2nd derivative of y.
   */
void
spline(double x[], double y[], int n, double yp1, double ypn, double y2[])
{
  int i, k;
  double p, qn, sig, un, *u;

#ifdef USE_ALLOCA
  u = (double*)alloca( (size_t) ((n-1)*sizeof(double)) );
#else
  u = vectorc(n-1);
#endif
  if (yp1 > 0.99e30)
    y2[0]=u[0]=0.0;
  else
    {
      y2[0] = -0.5;
      u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
    }
  for (i=1; i<=n-2; i++)
    {
      sig = (x[i]-x[i-1])/(x[i+1]-x[i-1]);
      p = sig*y2[i-1]+2.0;
      y2[i] = (sig-1.0)/p;
      u[i] = (y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
      u[i] = (6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
    }
  if (ypn > 0.99e30)
    qn=un=0.0;
  else
    {
      qn = 0.5;
      un = (3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
    }
  y2[n-1] = (un-qn*u[n-2])/(qn*y2[n-2]+1.0);
  for (k=n-2; k>=0; k--)
    y2[k] = y2[k]*y2[k+1]+u[k];
#ifndef USE_ALLOCA
  free_vectorc(u);
#endif
}
