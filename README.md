This is a set of functions for generating derivatives of the Cartesian
distance function $R=\sqrt{x^{2}+y^{2}+z^{2}}$. The main application
is intended to be numerical methods which use expansions of Green's
functions such as in the Laplace equation $G=1/4\pi R$.

# Prerequisites

The functions are deliberately kept simple so the files should compile
using any standard C compiler. There are two files, the main code
`radial.c` and a header `radial.h`. There is a third file `test.c`
which runs test code to check the installation. 

# Installation

The code is deliberately simple and intended to be integrated into
your own source code. To install the test program,

`make`

and to run it

`./radial-test`

# Usage

The main function is `radial_derivatives` which is called as:

`radial_derivatives(N, x, y, z, D) ;`

where $N$ is the highest order of derivative to be evaluated,
$(x,y,z)$ is the vector whose length $r$ is to be differentiated. On
exit, $D$ contains the scaled derivatives
$(1/\ell!m!n!)\partial^{\ell+m+n}/\partial x^{\ell}\partial
y^{m}\partial z^{n}$, i.e. the coefficients of the Taylor series for
$R$ expanded about $(x,y,z)$.

The location of derivatives in the output array are found using the
macros `radial_offset` and `radial_index`. Derivatives of total order
$N$ are located starting at index `radial_offset(N)` and derivatives
of order $(\ell,m,n)$ are found within that set of derivatives at
`radial_index(l,m,n)`. Thus the derivative
$(1/\ell!m!n!)\partial^{\ell+m+n}/\partial x^{\ell}\partial
y^{m}\partial z^{n}$ is located at

`D[i], i = radial_offset(l+m+n) + radial_index(l,m,n)`

The array $D$ must have at least `radial_offset(N+1)` elements. No
other workspace is required. 

