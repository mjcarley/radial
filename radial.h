#ifndef __RADIAL_H_INCLUDED__
#define __RADIAL_H_INCLUDED__

#define radial_offset(_n)      ((_n)*((_n)+1)*((_n)+2)/6)
#define radial_index(_l,_m,_n) ((_l)*(2*((_l)+(_m)+(_n))+3-(_l))/2+(_m))

int radial_derivatives(int N, double x, double y, double z, double *D) ;

#endif /*__RADIAL_H_INCLUDED__*/
