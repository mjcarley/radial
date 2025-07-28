#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#include "radial.h"

static double derivative_check(int N,
			       double x, double y, double z,
			       int d[], double ee, double *D,
			       double *D1, double *D2)
  
{
  int l, m, n, idx, jdx, off, q ;
  double dr, sc, err ;

  radial_derivatives(N-1, x+d[0]*ee/2, y+d[1]*ee/2, z+d[2]*ee/2, D1) ;
  radial_derivatives(N-1, x-d[0]*ee/2, y-d[1]*ee/2, z-d[2]*ee/2, D2) ;

  err = 0 ;
  for ( q = 0 ; q <= N-1 ; q ++ ) {
    off = radial_offset(q) ;
    for ( l = 0 ; l <= q ; l ++ ) {
      for ( m = 0 ; m <= q-l ; m ++ ) {
	n = q - l - m ;
	idx = radial_index(l,m,n) ;
	dr = (D1[off+idx] - D2[off+idx])/ee ;
	sc = tgamma(l+1)*tgamma(m+1)*tgamma(n+1)/
	  (tgamma(l+d[0]+1)*tgamma(m+d[1]+1)*tgamma(n+d[2]+1)) ;
	fprintf(stdout, "%d %d %d %lg ", l, m, n, dr*sc) ;
	jdx = radial_offset(q+1) + radial_index(l+d[0],m+d[1],n+d[2]) ;
	fprintf(stdout, "%lg\n", D[jdx]) ;
	if ( err < fabs(D[jdx] - dr*sc) ) err = fabs(D[jdx] - dr*sc) ;
      }
    }
  }
  
  return err ;
}

static void test_derivatives(int N, double x, double y, double z,
			     double *D, double *D1, double *D2, double ee)

{
  double err = 0 ;
  int d[3] ;

  fprintf(stderr, "derivatives w.r.t. x\n") ;
  d[0] = 1 ; d[1] = 0 ; d[2] = 0 ;
  err = derivative_check(N, x, y, z, d, ee, D, D1, D2) ;
  fprintf(stderr, "maximum error = %lg\n\n", err) ;

  fprintf(stderr, "derivatives w.r.t. y\n") ;
  d[0] = 0 ; d[1] = 1 ; d[2] = 0 ;
  err = derivative_check(N, x, y, z, d, ee, D, D1, D2) ;
  fprintf(stderr, "maximum error = %lg\n\n", err) ;

  fprintf(stderr, "derivatives w.r.t. z\n") ;
  d[0] = 0 ; d[1] = 0 ; d[2] = 1 ;
  err = derivative_check(N, x, y, z, d, ee, D, D1, D2) ;
  fprintf(stderr, "maximum error = %lg\n\n", err) ;

  return ;
}

static void test_differentiation(int N, double x, double y, double z,
				 int l, int m, int n,
				 double *D, double *D1, double *D2, double ee)

{
  double err = 0, de, dc, dx, dy, dz ;

  dx = 0.1 ; dy = -0.03 ; dz = 0.2 ;

  fprintf(stderr, "evaluation of derivatives\n") ;
  fprintf(stderr, "=========================\n") ;
  fprintf(stderr, "(dx,dy,dz) = (%lg,%lg,%lg)\n", dx, dy, dz) ;
  fprintf(stderr, "(l,m,n) = (%d,%d,%d)\n", l, m, n) ;
  
  de = radial_evaluate_derivative(N, D, dx, dy, dz, l, m, n) ;
  fprintf(stderr, "analytical derivative: %lg\n", de) ;

  if ( l > 0 ) {
    dc = (radial_evaluate_derivative(N, D, dx+ee/2, dy, dz, l-1, m, n) -
	  radial_evaluate_derivative(N, D, dx-ee/2, dy, dz, l-1, m, n))/ee ;

    fprintf(stderr, "numerical derivative (x): %lg\n", dc) ;
    if ( err < fabs(dc - de) ) err = fabs(dc - de) ;
  }

  if ( m > 0 ) {
    dc = (radial_evaluate_derivative(N, D, dx, dy+ee/2, dz, l, m-1, n) -
	  radial_evaluate_derivative(N, D, dx, dy-ee/2, dz, l, m-1, n))/ee ;

    fprintf(stderr, "numerical derivative (y): %lg\n", dc) ;
    if ( err < fabs(dc - de) ) err = fabs(dc - de) ;
  }
  
  if ( n > 0 ) {
    dc = (radial_evaluate_derivative(N, D, dx, dy, dz+ee/2, l, m, n-1) -
	  radial_evaluate_derivative(N, D, dx, dy, dz-ee/2, l, m, n-1))/ee ;

    fprintf(stderr, "numerical derivative (z): %lg\n", dc) ;
    if ( err < fabs(dc - de) ) err = fabs(dc - de) ;
  }

  fprintf(stderr, "maximum error: %lg\n", err) ;

  return ;
}

int main(int argc, char **argv)

{
  double *D, *D1, *D2, ee, x, y, z, dx, dy, dz, err, R, Rc ;
  int N, l, m, n, test ;
  char ch ;
  
  N = 40 ;

  x = 0.3 ; y = 1.9 ; z = -0.5 ;
  l = m = n = 0 ;
  ee = 1e-5 ;

  test = 2 ;
  
  while ( ( ch = getopt(argc, argv, "l:m:n:N:x:y:z:") ) != EOF ) {
    switch ( ch ) {
    default: fprintf(stderr, "%s: unrecognised option \"-%c\"\n",
		     argv[0], ch) ;
      return 1 ;
      break ;
    case 'l': l = atoi(optarg) ; break ;
    case 'm': m = atoi(optarg) ; break ;
    case 'n': n = atoi(optarg) ; break ;
    case 'N': N = atoi(optarg) ; break ;
    case 'x': x = atof(optarg) ; break ;
    case 'y': y = atof(optarg) ; break ;
    case 'z': z = atof(optarg) ; break ;
    }
  }
  
  D  = (double *)malloc(radial_offset(N+1)*sizeof(double)) ;
  D1 = (double *)malloc(radial_offset(N+1)*sizeof(double)) ;
  D2 = (double *)malloc(radial_offset(N+1)*sizeof(double)) ;

  radial_derivatives(N, x, y, z, D) ;

  fprintf(stderr, "maximum order: N = %d\n", N) ;
  fprintf(stderr, "evaluation point: (x,y,z) = (%lg,%lg,%lg)\n", x, y, z) ;
  fprintf(stderr, "radius: r = %lg\n", D[0]) ;

  if ( test == 0 ) {
    test_derivatives(N, x, y, z, D, D1, D2, ee) ;

    return 0 ;
  }
  
  if ( test == 1 ) {
    /*evaluation of expansion for distance*/
    err = 0.0 ;
    
    for ( dx = -0.5 ; dx <= 0.5 ; dx += 0.5/32 ) {
      for ( dy = -0.5 ; dy <= 0.5 ; dy += 0.5/32 ) {
      for ( dz = -0.5 ; dz <= 0.5 ; dz += 0.5/32 ) {
	R = radial_evaluate(N, D, dx, dy, dz) ;
	Rc = sqrt((x+dx)*(x+dx) + (y+dy)*(y+dy) + (z+dz)*(z+dz)) ;
	if ( fabs(Rc-R) > err ) err = fabs(Rc-R) ;
      }
      }
    }

    fprintf(stderr, "maximum error in distance: %lg\n", err) ;

    return 0 ;
  }

  if ( test == 2 ) {
    test_differentiation(N, x, y, z, l, m, n, D, D1, D2, ee) ;

    return 0 ;
  }  
  
  return 0 ;
}
