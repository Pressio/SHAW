
# Scope
C++ code for simulating elastic seismic shear waves in an axisymmetric domain.
The code uses Kokkos, and provides capabilities to run both 
the full order model (FOM) as well as the Galerkin reduced order model (ROM). 

# Building

## Prerequisites

You need to have BLAS installed. 
You need to have: 
```bash
export CC=<path-to-your-C-compiler>
export CXX=<path-to-your-CXX-compiler>
```

## Steps

- Build/install kokkos and kokkos-kernels 
- Use [this file](./do_build.sh) as follows: 
```bash
./do_build.sh \
 -working-dir=<where-you-want-to-build> \
 -kokkos-pfx=<path-to-your-kokkos-installation> \
 -kokkos-ker-pfx=<path-to-your-kokkos-kernels-installation> \
 --omp=yes
```



