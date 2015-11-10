#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#ifdef __unix__
#include <unistd.h>
#endif
#include "mynrutil.h"

#define NR_END 1
#define FREE_ARG void*


/* Library copies of inline functions for linkage with not optimized
   code. */
double
sqr(double a)
{
  return (a == 0.0 ? 0.0 : a*a);
}

float
fsqr(float a)
{
  return (a == 0.0 ? 0.0 : a*a);
}

double
dmax(double a, double b)
{
  return ( a>b ? a : b );
}

double
dmin(double a, double b)
{
  return ( a<b ? a : b );
}

long
lmax(long a, long b)
{
  return ( a>b ? a : b );
}

long
lmin(long a, long b)
{
  return ( a<b ? a : b );
}

int
imax(int a, int b)
{
  return ( a>b ? a : b );
}

int
imin(int a, int b)
{
  return ( a<b ? a : b );
}

double
sign(double a, double b)
{
  return ( b>=0.0 ? fabs(a) : -fabs(a));
}

/* Free a double vector allocated with vectorc(). */
#ifdef __GNUC__
void
free_vectorc(double *v)
{
  free((void*) v);
}
#endif


/* Numerical Recipes standard error handler. */
void
nrerror(char error_text[])
{
  fprintf(stderr, "Numerical Recipes run-time error...\n");
  fprintf(stderr, "%s\n", error_text);
  fprintf(stderr, "...now exiting to system...\n");
  exit(EXIT_FAILURE);
}

/* Sets time between attempts to malloc.  0 means exit on failure. */
unsigned int malloc_error_timeout = 0;

/* Report memory allocation error and exit if needed */
static void
malloc_error(char error_text[])
{
    fprintf(stderr, "Numerical Recipes run-time memory error...\n");
    fprintf(stderr, "%s\n", error_text);
    if ( malloc_error_timeout ) {
	fprintf( stderr,
		 "Waiting %d sec for memory...\n", malloc_error_timeout );
	sleep(malloc_error_timeout);
    } else {
	nrerror("memory allocation error");
    }
}


/* Allocate a double vector with subscript range v[0..n-1]. */
double *
vectorc(long n)
{
  double *v;
  v = (double*) malloc((size_t) (n*sizeof(double)));
  while (!v) {
      malloc_error("allocation failure in vectorc()");
      v = (double*) malloc((size_t) (n*sizeof(double)));
  }
  return v;
}


/* Allocate a double vector with subscript range v[nl..nh]. */
double *
vector(long nl, long nh)
{
  double *v;

  v = (double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
  while (!v) {
      malloc_error("allocation failure in vector()");
      v = (double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
  }
  return v-nl+NR_END;
}

/* Allocate an int vector with subscript range v[nl..nh]. */
int *
ivector(long nl, long nh)
{
  int *v;

  v = (int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
  while (!v) {
      malloc_error("allocation failure in ivector()");
      v = (int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
  }
  return v-nl+NR_END;
}

/* Allocate an unsigned char vector with subscript range v[nl..nh]. */
unsigned char *
cvector(long nl, long nh)
{
  unsigned char *v;

  v = (unsigned char *) malloc((size_t) 
			       ((nh-nl+1+NR_END)*sizeof(unsigned char)));
  while (!v) {
      malloc_error("allocation failure in cvector()");
      v = (unsigned char *) malloc((size_t) 
				   ((nh-nl+1+NR_END)*sizeof(unsigned char)));
  }
  return v-nl+NR_END;
}

/* Allocate an unsigned long vector with subscript range v[nl..nh]. */
unsigned long *
lvector(long nl, long nh)
{
  unsigned long *v;

  v = (unsigned long *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(long)));
  while (!v) {
      malloc_error("allocation failure in lvector()");
      v = (unsigned long *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(long)));
  }
  return v-nl+NR_END;
}

/* Allocate a float vector with subscript range v[nl..nh]. */
float *
fdvector(long nl, long nh)
{
  float *v;

  v = (float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
  while (!v) {
      malloc_error("allocation failure in fvector()");
      v = (float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
  }
  return v-nl+NR_END;
}

/* Allocate a double matrix with subscript range
   m[nrl..nrh][ncl..nch]. */
double **
matrix(long nrl, long nrh, long ncl, long nch)
{
  long i, nrow=nrh-nrl+1, ncol=nch-ncl+1;
  double **m;

  /* allocate pointers to rows */
  m = (double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
  while (!m) {
      malloc_error("allocation failure 1 in matrix()");
      m = (double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
  }
  m += NR_END;
  m -= nrl;

  /* allocate rows and set pointers to them */
  m[nrl] = (double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
  while (!m[nrl]) {
      malloc_error("allocation failure 2 in matrix()");
      m[nrl] = (double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
  }
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for(i=nrl+1; i<=nrh; i++)
    m[i] = m[i-1]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

/* Allocate a float matrix with subscript range
   m[nrl..nrh][ncl..nch]. */
float **
fmatrix(long nrl, long nrh, long ncl, long nch)
{
  long i, nrow=nrh-nrl+1, ncol=nch-ncl+1;
  float **m;

  /* allocate pointers to rows */
  m = (float **) malloc((size_t)((nrow+NR_END)*sizeof(float*)));
  while (!m) {
      malloc_error("allocation failure 1 in fmatrix()");
      m = (float **) malloc((size_t)((nrow+NR_END)*sizeof(float*)));
  }
  m += NR_END;
  m -= nrl;

  /* allocate rows and set pointers to them */
  m[nrl] = (float *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float)));
  while (!m[nrl]) {
      malloc_error("allocation failure 2 in fmatrix()");
      m[nrl] = (float *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float)));
  }
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for(i=nrl+1; i<=nrh; i++)
    m[i] = m[i-1]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

/* Allocate a int matrix with subscript range
   m[nrl..nrh][ncl..nch]. */
int **
imatrix(long nrl, long nrh, long ncl, long nch)
{
  long i, nrow=nrh-nrl+1, ncol=nch-ncl+1;
  int **m;

  /* allocate pointers to rows */
  m = (int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
  while (!m) {
      malloc_error("allocation failure 1 in imatrix()");
      m = (int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
  }
  m += NR_END;
  m -= nrl;

  /* allocate rows and set pointers to them */
  m[nrl] = (int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
  while (!m[nrl]) {
      malloc_error("allocation failure 2 in imatrix()");
      m[nrl] = (int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
  }
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for(i=nrl+1; i<=nrh; i++)
    m[i] = m[i-1]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

/* Point a submatrix [newrl..][newcl..] to
   a[oldrl..oldrh][oldcl..oldch] */
double **
submatrix(double **a, long oldrl, long oldrh, long oldcl, long oldch,
	  long newrl, long newcl)
{
  long i, j, nrow=oldrh-oldrl+1, ncol=oldcl-newcl;
  double **m;

  /* allocate array of pointers to rows */
  m = (double **) malloc((size_t) ((nrow+NR_END)*sizeof(double*)));
  while (!m) {
      malloc_error("allocation failure in submatrix()");
      m = (double **) malloc((size_t) ((nrow+NR_END)*sizeof(double*)));
  }
  m += NR_END;
  m -= newrl;

  /* set pointers to rows */
  for(i=oldrl,j=newrl; i<=oldrh; i++,j++)
    m[j] = a[i]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

/* Allocate a double matrix m[nrl..nrh][ncl..nch] that points to the
   matrix declared in the standard C manner as a[nrow][ncol], where
   nrow=nrh-nrl+1 and ncol=nch-ncl+1.  The routine should be called
   with the address &a[0][0] as the first argument. */
double **
convert_matrix(double *a, long nrl, long nrh, long ncl, long nch)
{
  long i, j, nrow=nrh-nrl+1, ncol=nch-ncl+1;
  double **m;

  /* allocate pointers to rows */
  m = (double **) malloc((size_t) ((nrow+NR_END)*sizeof(double*)));
  while (!m) {
      malloc_error("allocation failure in convert_matrix()");
      m = (double **) malloc((size_t) ((nrow+NR_END)*sizeof(double*)));
  }
  m += NR_END;
  m -= nrl;

  /* set pointers to rows */
  m[nrl] = a-ncl;
  for(i=1,j=nrl+1; i<nrow; i++,j++)
    m[j] = m[j-1]+ncol;
  /* return pointer to array of pointers to rows */
  return m;
}

/* Allocate a double 3tensor with range
   t[nrl..nrh][ncl..nch][ndl..ndh]. */
double ***
d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
{
  long i, j, nrow=nrh-nrl+1, ncol=nch-ncl+1, ndep=ndh-ndl+1;
  double ***t;

  /* allocate pointers to pointers to rows */
  t = (double ***) malloc((size_t)((nrow+NR_END)*sizeof(double**)));
  while (!t) {
      malloc_error("allocation failure 1 in d3tensor()");
      t = (double ***) malloc((size_t)((nrow+NR_END)*sizeof(double**)));
  }
  t += NR_END;
  t -= nrl;

  /* allocate pointers to rows and set pointers to them */
  t[nrl] = (double **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double*)));
  while (!t[nrl]) {
      malloc_error("allocation failure 2 in d3tensor()");
      t[nrl] = (double **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double*)));
  }
  t[nrl] += NR_END;
  t[nrl] -= ncl;

  /* allocate rows and set pointers to them */
  t[nrl][ncl] = ((double *)
		 malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(double))));
  while (!t[nrl][ncl]) {
      malloc_error("allocation failure 3 in d3tensor()");
      t[nrl][ncl] = ((double *)
		     malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(double))));
  }
  t[nrl][ncl] += NR_END;
  t[nrl][ncl] -= ndl;

  for(j=ncl+1; j<=nch; j++)
    t[nrl][j] = t[nrl][j-1]+ndep;
  for(i=nrl+1; i<=nrh; i++)
    {
      t[i] = t[i-1]+ncol;
      t[i][ncl] = t[i-1][ncl]+ncol*ndep;
      for(j=ncl+1; j<=nch; j++)
	t[i][j]=t[i][j-1]+ndep;
    }

  /* return pointer to array of pointers to rows */
  return t;
}



/* Free a double vector allocated with vector(). */
void
free_vector(double *v, long nl, long nh)
{
  free((FREE_ARG) (v+nl-NR_END));
}

/* Free an int vector allocated with ivector(). */
void
free_ivector(int *v, long nl, long nh)
{
  free((FREE_ARG) (v+nl-NR_END));
}

/* Free an unsigned char vector allocated with cvector(). */
void
free_cvector(unsigned char *v, long nl, long nh)
{
  free((FREE_ARG) (v+nl-NR_END));
}

/* Free an unsigned long vector allocated with lvector(). */
void
free_lvector(unsigned long *v, long nl, long nh)
{
  free((FREE_ARG) (v+nl-NR_END));
}

/* Free a float vector allocated with fvector(). */
void
free_fvector(float *v, long nl, long nh)
{
  free((FREE_ARG) (v+nl-NR_END));
}

/* Free a double matrix allocated by matrix(). */
void
free_matrix(double **m, long nrl, long nrh, long ncl, long nch)
{
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}

/* Free a float matrix allocated by fmatrix(). */
void
free_dmatrix(float **m, long nrl, long nrh, long ncl, long nch)
{
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}

/* Free an int matrix allocated by imatrix(). */
void
free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
{
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}

/* Free a submatrix allocated by submatrix(). */
void
free_submatrix(double **b, long nrl, long nrh, long ncl, long nch)
{
  free((FREE_ARG) (b+nrl-NR_END));
}

/* Free a matrix allocated by convert_matrix(). */
void
free_convert_matrix(double **b, long nrl, long nrh, long ncl, long nch)
{
  free((FREE_ARG) (b+nrl-NR_END));
}

/* Free a double d3tensor allocated by d3tensor(). */
void
free_d3tensor(double ***t, long nrl, long nrh, long ncl, long nch,
	      long ndl, long ndh)
{
  free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
  free((FREE_ARG) (t[nrl]+ncl-NR_END));
  free((FREE_ARG) (t+nrl-NR_END));
}
