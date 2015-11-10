#ifndef _NR_UTILS_H_
#define _NR_UTILS_H_

#include <math.h>
#include <stdlib.h>

#ifndef __GNUC__
#define __inline__
#define __const__
#endif


#ifdef __cplusplus
#  define __BEGIN_DECLS extern "C" {
#  define __END_DECLS }
#else
#  define __BEGIN_DECLS
#  define __END_DECLS
#endif

#ifdef __GNUC__
extern __inline__ double
sqr(double a)
{
  return (a == 0.0 ? 0.0 : a*a);
}
#define SQR(a) (sqr(a))

extern __inline__ float
fsqr(float a)
{
  return (a == 0.0 ? 0.0 : a*a);
}
#define FSQR(a) (rsqr(a))

extern __inline__ double
dmax(double a, double b)
{
  return ( a>b ? a : b );
}
#define DMAX(a,b) (dmax((a),(b)))

extern __inline__ double
dmin(double a, double b)
{
  return ( a<b ? a : b );
}
#define DMIN(a,b) (dmin((a),(b)))

extern __inline__ double
fmax(double a, double b)
{
  return ( a>b ? a : b );
}
#define FMAX(a,b) (fmax((a),(b)))

extern __inline__ double
fmin(double a, double b)
{
  return ( a<b ? a : b );
}
#define FMIN(a,b) (fmin((a),(b)))
extern __inline__ long
lmax(long a, long b)
{
  return ( a>b ? a : b );
}
#define LMAX(a,b) (lmax((a),(b)))

extern __inline__ long
lmin(long a, long b)
{
  return ( a<b ? a : b );
}
#define LMIN(a,b) (lmin((a),(b)))

extern __inline__ int
imax(int a, int b)
{
  return ( a>b ? a : b );
}
#define IMAX(a,b) (imax((a),(b)))

extern __inline__ int
imin(int a, int b)
{
  return ( a<b ? a : b );
}
#define IMIN(a,b) (imin((a),(b)))

extern __inline__ double
sign(double a, double b)
{
  return ( b>=0.0 ? fabs(a) : -fabs(a));
}
#define SIGN(a,b) (sign((a),(b)))
#else /* !__GNUC__ */
extern double sqr(double a);
#define SQR(a) (sqr(a))
extern float fsqr(float a);
#define FSQR(a) (rsqr(a))
extern double dmax(double a, double b);
#define DMAX(a,b) (dmax((a),(b)))
extern double dmin(double a, double b);
#define DMIN(a,b) (dmin((a),(b)))
#if 0
extern __inline__ double fmax(double a, double b);
#define FMAX(a,b) (fmax((a),(b)))
extern __inline__ double fmin(double a, double b);
#define FMIN(a,b) (fmin((a),(b)))
#endif
extern __inline__ long lmax(long a, long b);
#define LMAX(a,b) (lmax((a),(b)))
extern __inline__ long lmin(long a, long b);
#define LMIN(a,b) (lmin((a),(b)))
extern __inline__ int imax(int a, int b);
#define IMAX(a,b) (imax((a),(b)))
extern __inline__ int imin(int a, int b);
#define IMIN(a,b) (imin((a),(b)))
extern __inline__ double sign(double a, double b);
#define SIGN(a,b) (sign((a),(b)))
#endif /* !__GNUC__ */

__BEGIN_DECLS

extern unsigned int malloc_error_timeout;

void nrerror(char error_text[]);
double *vector(long nl, long nh);
int *ivector(long nl, long nh);
unsigned char *cvector(long nl, long nh);
unsigned long *lvector(long nl, long nh);
float *fvector(long nl, long nh);
double **matrix(long nrl, long nrh, long ncl, long nch);
float **fmatrix(long nrl, long nrh, long ncl, long nch);
int **imatrix(long nrl, long nrh, long ncl, long nch);
double **submatrix(double **a, long oldrl, long oldrh, long oldcl, long oldch,
		   long newrl, long newcl);
double **convert_matrix(double *a, long nrl, long nrh, long ncl, long nch);
double ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);

double *vectorc(long n);
/* Free a double vector allocated with vectorc(). */
#ifdef __GNUC__
extern __inline__ void
free_vectorc(double *v)
{
  free((void*) v);
}
#else
#define free_vectorc(v) (free((void*)(v)))
#endif

void free_vector(double *v, long nl, long nh);
void free_ivector(int *v, long nl, long nh);
void free_cvector(unsigned char *v, long nl, long nh);
void free_lvector(unsigned long *v, long nl, long nh);
void free_fvector(float *v, long nl, long nh);
void free_matrix(double **m, long nrl, long nrh, long ncl, long nch);
void free_fmatrix(float **m, long nrl, long nrh, long ncl, long nch);
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);
void free_submatrix(double **b, long nrl, long nrh, long ncl, long nch);
void free_convert_matrix(double **b, long nrl, long nrh, long ncl, long nch);
void free_d3tensor(double ***t, long nrl, long nrh, long ncl, long nch,
	long ndl, long ndh);

__END_DECLS

#endif /* _NR_UTILS_H_ */
