v2DM-DOCI
=========
This program implements a v2DM solver for the DOCI N-representability constraints.
Both the boundary point method and the potential reduction method have been
implemented. Two orbital optimizer are available:
- A simulated annealing code (in the extern directory)
- Jacobi rotations based local minimalisation

Build
-----
To build this, you need a C++11 compiler (gcc 4.8 or newer, clang 3.3 or newer),
the HDF5 libraries and the blas and lapack libraries. The Makefile is quite 
simple, adjust the compilers and header/libraries as needed for your system.

Input
-----
The program needs symmetry adapted electron integrals from [PSI4](https://github.com/psi4/psi4public). 
You can extract them using a plugin found in the extern/mointegrals.cpp file. 

Documentation
-------------
Everything is documented with doxygen. Run `make doc` to build the HTML docs. 
Or read the docs [online](http://wpoely86.github.io/doci_sdp-atom/). If something
is unclear, do not hesitate to contact me.

License
-------
The code is available under the [GPLv3](https://www.gnu.org/licenses/gpl-3.0.txt) license.
