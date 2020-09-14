
# Scope
This repository contains code for simulating elastic seismic shear waves in an axisymmetric domain.
The implementation uses Kokkos and Kokkos-kernels to enable performance portability,
and we provide capabilities to run both the full order model (FOM)
as well as the Galerkin reduced order model (ROM).

Specifically, this code has been developed mainly for a paper submitted on 11th Sept. 2020 to SISC which is under review.
(We will add a link as soon as we get the paper back or we post it on arXiv).

# Prerequisites
To build and use the code, you need to have CMake>=3.11.0,
C, C++, Fortran compilers, as well as Python>3.6 with at least
the following packages: matplotlib, numpy, scipy, yaml.
The code has been tested on MaxOS Catalina 10.15.5 and RedHat 7.0
with GCC-8.3.1 and GCC-8.4.0.

# Building
A [step-by-step guide](./docs/build.md) on how to build/install the
TPLs needed and the elastic shear wave code.

# Step-by-step for running a full-order model simulation
Follow [this guide](./docs/run_fom.md) for an example to run a FOM run.

# Step-by-step for a rank-1 Galerkin ROM simulation
Follow [this guide](./docs/run_rom.md) for an example to run a rank-1 Galerkin run.

# Step-by-step for a rank-2 Galerkin ROM simulation
Coming soon.
<!--
Follow [this guide](./docs/run_rom2.md) for an example to run a rank-2 Galerkin
to efficiently generate samples of the forcing period. -->

# Content
- [bash script driving the build](./do_build.sh)
- [C++ source code and tests](./cpp)
- [meshing](./meshing)
- [Python scripts for processing and workflows](./python_scripts)
