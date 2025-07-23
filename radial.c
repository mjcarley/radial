#include <stdio.h>
#include <math.h>

#include "radial.h"

int radial_derivatives(int N, double x, double y, double z, double *D)

{
  int l, m, n, q, idx, jdx, kdx, off, i1, i2, s, t, u ;
  double r ;

  off = radial_offset(N+1) ;
  for ( n = 0 ; n < off ; n ++ ) D[n] = 0.0 ;
  
  r = sqrt(x*x + y*y + z*z) ;
  off = radial_offset(0) ;
  idx = radial_index(0,0,0) ;
  D[off+idx] = r ;

  if ( N == 0 ) return 0 ;

  /*generate single variable derivatives*/
  off = radial_offset(1) ;
  idx = radial_index(1,0,0) ; D[off+idx] = x/r ;  
  idx = radial_index(0,1,0) ; D[off+idx] = y/r ;  
  idx = radial_index(0,0,1) ; D[off+idx] = z/r ;  

  if ( N == 1 ) return 0 ;

  /*easier to include mixed second-order derivatives here in case N==2*/
  off = radial_offset(2) ;
  idx = radial_index(2,0,0) ; D[off+idx] = (1.0 - x*x/r/r)/r/2.0 ;
  idx = radial_index(1,1,0) ; D[off+idx] = -x*y/r/r/r ;
  idx = radial_index(0,2,0) ; D[off+idx] = (1.0 - y*y/r/r)/r/2.0 ;  
  idx = radial_index(0,1,1) ; D[off+idx] = -y*z/r/r/r ;
  idx = radial_index(0,0,2) ; D[off+idx] = (1.0 - z*z/r/r)/r/2.0 ;  
  idx = radial_index(1,0,1) ; D[off+idx] = -z*x/r/r/r ;
  
  if ( N == 2 ) return 0 ;

  /*unmixed derivatives*/
  for ( n = 3 ; n <= N ; n ++ ) {
    idx = radial_offset(n) + radial_index(n, 0, 0) ;
    jdx = radial_offset(n) + radial_index(0, n, 0) ;
    kdx = radial_offset(n) + radial_index(0, 0, n) ;

    for ( l = 1 ; l <= n-1 ; l ++ ) {
      i1 = radial_offset(  l) + radial_index(  l,   0,   0) ;
      i2 = radial_offset(n-l) + radial_index(n-l,   0,   0) ;
      D[idx] += D[i1]*D[i2] ;
      i1 = radial_offset(  l) + radial_index(  0,   l,   0) ;
      i2 = radial_offset(n-l) + radial_index(  0, n-l,   0) ;
      D[jdx] += D[i1]*D[i2] ;
      i1 = radial_offset(  l) + radial_index(  0,   0,   l) ;
      i2 = radial_offset(n-l) + radial_index(  0,   0, n-l) ;
      D[kdx] += D[i1]*D[i2] ;
    }
    D[idx] /= -2*r ; D[jdx] /= -2*r ; D[kdx] /= -2*r ;
  }
  
  for ( q = 3 ; q <= N ; q ++ ) {
    /*two-variable derivatives*/
    l = 0 ;
    for ( m = 1 ; m <= q-1 ; m ++ ) {
      n = q - m - l ;
      idx = radial_offset(l+m+n) + radial_index(l, m, n) ;
      for ( t = 0 ; t <= n-1 ; t ++ ) {
	i1 = radial_offset(n-t) + radial_index(0,0,n-t) ;
	i2 = radial_offset(m+t) + radial_index(0, m, t) ;
	D[idx] += 2.0*D[i1]*D[i2] ;
      /* } */
      /* for ( t = 0 ; t <= n ; t ++ ) { */
	for ( s = 1 ; s <= m-1 ; s ++ ) {
	  i1 = radial_offset(s+n-t) + radial_index(0, s, n-t) ;
	  i2 = radial_offset(m-s+t) + radial_index(0, m-s, t) ;
	  D[idx] += D[i1]*D[i2] ;
	}
      }
      t = n ;
      for ( s = 1 ; s <= m-1 ; s ++ ) {
	i1 = radial_offset(s+n-t) + radial_index(0, s, n-t) ;
	i2 = radial_offset(m-s+t) + radial_index(0, m-s, t) ;
	D[idx] += D[i1]*D[i2] ;
      }
      
      D[idx] /= -2*r ;
    }

    m = 0 ;
    for ( n = 1 ; n <= q - 1 ; n ++ ) {
      l = q - m - n ;
      idx = radial_offset(l+m+n) + radial_index(l, m, n) ;
      for ( t = 0 ; t <= l-1 ; t ++ ) {
	i1 = radial_offset(l-t) + radial_index(l-t, 0, 0) ;
	i2 = radial_offset(n+t) + radial_index(t, 0, n) ;
	D[idx] += 2.0*D[i1]*D[i2] ;
      /* } */
      /* for ( t = 0 ; t <= l ; t ++ ) { */
	for ( s = 1 ; s <= n-1 ; s ++ ) {
	  i1 = radial_offset(l-t+s) + radial_index(l-t, 0, s) ;
	  i2 = radial_offset(n-s+t) + radial_index(t, 0, n-s) ;
	  D[idx] += D[i1]*D[i2] ;
	}
      }
      t = l ;
      for ( s = 1 ; s <= n-1 ; s ++ ) {
	i1 = radial_offset(l-t+s) + radial_index(l-t, 0, s) ;
	i2 = radial_offset(n-s+t) + radial_index(t, 0, n-s) ;
	D[idx] += D[i1]*D[i2] ;
      }
      D[idx] /= -2*r ;
    }

    for ( l = 1 ; l <= q ; l ++ ) {
      for ( m = 1 ; m <= q-l ; m ++ ) {
	n = q - m - l ;
	idx = radial_offset(l+m+n) + radial_index(l, m, n) ;
      
	for ( u = 0 ; u <= n ; u ++ ) {
	  for ( t = 0 ; t <= m-1 ; t ++ ) {
	    for ( s = 1 ; s <= l-1 ; s ++ ) {
	      i1 = radial_offset(s+m-t+u) + radial_index(s, m-t, u) ;
	      i2 = radial_offset(l-s+t+n-u) + radial_index(l-s, t, n-u) ;
	      D[idx] += D[i1]*D[i2] ;
	    }
	    /* } */
	    /* for ( t = 0 ; t <= m-1 ; t ++ ) { */
	    /*this is the s==0 case from the previous loop, but not the 
	      factor of 2*/
	    i1 = radial_offset(m-t+u) + radial_index(0, m-t, u) ;
	    i2 = radial_offset(l+t+n-u) + radial_index(l, t, n-u) ;
	    D[idx] += 2.0*D[i1]*D[i2] ;
	  }
	  t = m ;
	  for ( s = 1 ; s <= l-1 ; s ++ ) {
	    i1 = radial_offset(s+m-t+u) + radial_index(s, m-t, u) ;
	    i2 = radial_offset(l-s+t+n-u) + radial_index(l-s, t, n-u) ;
	    D[idx] += D[i1]*D[i2] ;
	  }
	  /* for ( u = 0 ; u <= n-1 ; u ++ ) { */
	  /* i1 = radial_offset(n-u) + radial_index(0,0,n-u) ; */
	  /* i2 = radial_offset(l+m+u) + radial_index(l,m,u) ; */
	  /* D[idx] += 2*D[i1]*D[i2] ; */
	  /* } */
	}
	/* u = n ; */
	/* 	i1 = radial_offset(n-u) + radial_index(0,0,n-u) ; */
	/* 	i2 = radial_offset(l+m+u) + radial_index(l,m,u) ; */
	/* 	D[idx] -= 2*D[i1]*D[i2] ; */
	for ( u = 0 ; u <= n-1 ; u ++ ) {
	  i1 = radial_offset(n-u) + radial_index(0,0,n-u) ;
	  i2 = radial_offset(l+m+u) + radial_index(l,m,u) ;
	  D[idx] += 2*D[i1]*D[i2] ;
	}
	D[idx] /= -2*r ;
      }
    }
  }
	  
  return 0 ;
}
