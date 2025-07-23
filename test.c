#include <stdio.h>
#include <math.h>

#include "radial.h"

static double derivative_check(int N,
			       double x, double y, double z,
			       int d[], double ee, double *D)
  
{
  int l, m, n, idx, jdx, off, q ;
  double D1[8192], D2[8192], dr, sc, err ;

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

int main(int argc, char **argv)

{
  double D[8192],ee, x, y, z, err ;
  int N, d[3] ;

  N = 32 ;

  x = 0.3 ; y = 1.9 ; z = -0.5 ;
  ee = 1e-5 ;
  
  radial_derivatives(N, x, y, z, D) ;

  fprintf(stderr, "r = %lg\n", D[0]) ;

  fprintf(stderr, "derivatives w.r.t. x\n") ;
  d[0] = 1 ; d[1] = 0 ; d[2] = 0 ;
  err = derivative_check(N, x, y, z, d, ee, D) ;
  fprintf(stderr, "maximum error = %lg\n\n", err) ;
  fprintf(stderr, "derivatives w.r.t. y\n") ;
  d[0] = 0 ; d[1] = 1 ; d[2] = 0 ;
  err = derivative_check(N, x, y, z, d, ee, D) ;
  fprintf(stderr, "maximum error = %lg\n\n", err) ;

  fprintf(stderr, "derivatives w.r.t. z\n") ;
  d[0] = 0 ; d[1] = 0 ; d[2] = 1 ;
  err = derivative_check(N, x, y, z, d, ee, D) ;
  fprintf(stderr, "maximum error = %lg\n\n", err) ;

  return 0 ;
}
