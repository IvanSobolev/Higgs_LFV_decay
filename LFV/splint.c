#include "mynr.h"
#include "mynrutil.h"

/* Cubic spline interpolation.

   Input:
   xa --- vector of x values.
   ya --- vector of y values.
   y2a --- vector of 2nd derivative generated by spline.
   n --- number of data points (vector length).
   x --- x value.
   Output:
   *y --- y(x) interpolation.
   */
void
splint(double xa[], double ya[], double y2a[], int n, double x, double *y)
{
  int klo, khi, k;
  double h, b, a;

  klo = 0;
  khi = n-1;
  while (khi-klo > 1)
    {
      k = (khi+klo) >> 1;
      if (xa[k] > x)
	khi = k;
      else
	klo = k;
    }
  h = xa[khi]-xa[klo];
  if (h == 0.0)
    nrerror("Bad xa input to routine splint");
  a = (xa[khi]-x)/h;
  b = (x-xa[klo])/h;
  *y = a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}
