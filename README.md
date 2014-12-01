v2DM DOCI
=========
This program implements a v2DM solver for DOCI. Both the boundary point method and the
potential reduction method have been implemented. A simulated annealing code has been
added to do the geometric optimalisation (in the extern directory).

Build
-----
To build this, you need a C++11 compiler (gcc 4.8 or newer, clang 3.3 or newer),
the HDF5 libraries and the blas and lapack libraries. The
Makefile is quite simple, adjust the compilers and header/libraries as needed
for your system.

Input
-----
The program needs molecular integrals from [PSI4](https://github.com/psi4/psi4public). 
You can extract them using a plugin found in the extern/mointegrals.cpp file. 

Documentation
-------------
Everything is documented with doxygen. Run `make doc` to build the HTML docs.
