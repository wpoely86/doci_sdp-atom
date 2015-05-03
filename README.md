[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.17326.svg)](http://dx.doi.org/10.5281/zenodo.17326)

v2DM-DOCI
=========
This program implements a v2DM solver for the DOCI N-representability constraints.
Both the boundary point method and the potential reduction method have been
implemented. Two orbital optimizer are available:
- A simulated annealing code (in the extern directory)
- Jacobi rotations based local minimalisatie

Build
-----
To build this, you need a C++11 compiler (GCC 4.8 or newer, Clang 3.3 or newer),
the HDF5 libraries and the blas and lapack libraries. The Makefile is quite 
simple, adjust the compilers and header/libraries as needed for your system.

Input
-----
The program needs symmetry-adapted electron integrals from [PSI4](https://github.com/psi4/psi4public). 
You can extract them using a plugin found in the extern/mointegrals.cpp file. 

Documentation
-------------
Everything is documented with doxygen. Run `make doc` to build the HTML docs. 
Or read the docs [online](http://wpoely86.github.io/doci_sdp-atom/). If something
is unclear, do not hesitate to contact me.

Example
-------
In the directory `example`, all the input files for the He dimer (interatomic distance is 10 Bohr)
can be found:
- `He2-integrals.dat`: This is the input file for PSI4 to generate the one- and two-electron integrals
- `mo-ints-he2-cc-pvdz-010.00.h5`: the Hartree-Fock molecular orbitals
- `so-ints-he2-cc-pvdz-010.00.h5`: the integrals in the orthogonalized, symmetry-adapted basis
- `unitary-mo-he2-cc-pvdz-010.00.h5`: the orthogonal transformation from the orthogonalized, 
  symmetry-adapted basis to the Hartree-Fock molecular orbital.

You can run the v2DM-DOCI optimisation with: 
`../doci_bp -i so-ints-he2-cc-pvdz-010.00.h5 -u unitary-mo-he2-cc-pvdz-010.00.h -l`

The will generate:
- `optimale-uni.h5`: the optimal basis transformation
- `optimal-ham.h5`: the one- and two-electron integrals in the optimal basis
- `optimal-rdm.h5`: the optimal 2DM

In the file `output.txt`, the output of the optimisation is stored.

License
-------
The code is available under the [GPLv3](https://www.gnu.org/licenses/gpl-3.0.txt) license.
