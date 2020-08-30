
# Scope
C++ code for simulating elastic seismic shear waves in an axisymmetric domain.
The implementation uses Kokkos and Kokkos-kernels to enable performance portability,
and we provide capabilities to run both the full order model (FOM)
as well as the Galerkin reduced order model (ROM).

This code has been developed for: cite-article.
You can find more details in that paper.

# Content

- [bash script driving the build](./do_build.sh)
- [C++ source code and tests](./cpp)
- [meshing](./meshing)
- [Python scripts for processing and workflows](./python_scripts)

# Prerequisites
To build and use the code, you need to have CMake, serial C, C++
and Fortran compilers, as well as Python>3.6 with at least
the following packages: matplotlib, numpy, scipy, yaml.

# Building the code
A [step-by-step guide](./docs/build.md) on how to build/install all
TPLs needed and the elastic shear wave code.

# Creating and running a full-order model simulation
A [step-by-step guide](./docs/run.md) on how to setup a FOM run.
