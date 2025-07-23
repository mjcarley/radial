This is a set of functions for generating derivatives of the Cartesian
distance function $R=\sqrt{x^{2}+y^{2}+z^{2}}$. The main application
is intended to be numerical methods which use expansions of Green's
functions such as in the Laplace equation $G=1/4\pi R$.

# Prerequisites

The functions are deliberately kept simple so the files should compile
using any standard C compiler. There are two files, the main code
`radial.c` and a header `radial.h`. There is a third file `test.c`
which runs test code to check the installation. 

# Usage

The main function is `radial_derivatives`
