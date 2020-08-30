
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


# Building
Here is a step-by-step guide on how to build/install all
TPLs needed and the code.

## Step 1: set basic environment
To simplify the process, define the following env variables:
```bash
export ESWSRCDIR=<path-to-the-code-repository>
export CC=<path-to-your-C-compiler>
export CXX=<path-to-your-C++-compiler>
export FC=<path-to-your-Fortran-compiler>
export MYWORKDIR=<path-to-where-you-want-to-work-in>
```
and execute:
```bash
mkdir -p ${MYWORKDIR}
```

## Step 2: BLAS/LAPACK
Here we build and install BLAS and LAPACK.
Proceed as follows:
```bash
cd ${MYWORKDIR}
mkdir tpls && cd tpls
cp ${ESWSRCDIR}/bash/build_openblas.sh .
bash build_openblas.sh
```
and then define:
```bash
export BLAS_ROOT=${MYWORKDIR}/tpls/openblas/install
export BLASLIBNAME=openblas
export LAPACK_ROOT=${MYWORKDIR}/tpls/openblas/install
export LAPACKLIBNAME=openblas
```
this should build and install OpenBLAS such that
inside `$MYWORKDIR}/tpls/openblas/install/lib` you should see something as:
```bash
total 78328
drwxr-xr-x  3 fnrizzi  staff    96B Aug 30 09:40 cmake
lrwxr-xr-x  1 fnrizzi  staff    34B Aug 30 09:40 libopenblas.0.dylib
lrwxr-xr-x  1 fnrizzi  staff    30B Aug 30 09:40 libopenblas.a
lrwxr-xr-x  1 fnrizzi  staff    34B Aug 30 09:40 libopenblas.dylib
drwxr-xr-x  3 fnrizzi  staff    96B Aug 30 09:40 pkgconfig
```
**Note**: if you already have BLAS/LAPACK installed, you can skip
the build step above and directly set the needed env vars to
point to your BLAS/LAPACK installation.


## Step 3: Kokkos and Kokkos-kernels
Now that we have BLAS/LAPACK built, we build Kokkos and Kokkos-kernels.
Proceed as follows:
```bash
cd ${MYWORKDIR}/tpls
cp ${ESWSRCDIR}/bash/build_kokkos_and_kernels.sh .
export KOKKOSPFX=${MYWORKDIR}/tpls/kokkos/kokkos_install
export KOKKOSKERPFX=${MYWORKDIR}/tpls/kokkos/kokkos_kernels_install
bash build_kokkos_and_kernels.sh
```
**Remarks**:
* the above process builds Kokkos/Kokkos-kernels *without* any arch-specific
optimization, since this is meant to work on any system. However, if you want to
have arch-specific optimizations (and you should), you need to change the arch flag
passed to Kokkos (see inside `build_kokkos_and_kernels.sh`) and rebuild;
* we only enable the OpenMP backend;
* if you already have Kokkos/Kokkos-kernels installed, you can skip the build step
above and directly set the needed env vars to point to your installation.


## Step 4: Build the Elastic Shear Wave executables
Proceed as follows:
```bash
cd $ESWSRCDIR}
./do_build.sh \
 -working-dir=${MYWORKDIR} \
 -kokkos-pfx=${KOKKOSPFX} \
 -kokkos-ker-pfx=${KOKKOSKERPFX} \
 --omp=yes
```
this should generate inside `${MYWORKDIR}/build` the following executables:
```bash
-rw-r--r--   1 fnrizzi  staff    14K Aug 30 10:41 CMakeCache.txt
drwxr-xr-x  18 fnrizzi  staff   576B Aug 30 10:42 CMakeFiles
-rw-r--r--   1 fnrizzi  staff   333B Aug 30 10:41 CTestTestfile.cmake
-rw-r--r--   1 fnrizzi  staff    15K Aug 30 10:41 Makefile
-rw-r--r--   1 fnrizzi  staff   1.6K Aug 30 10:41 cmake_install.cmake
-rwxr-xr-x   1 fnrizzi  staff    38K Aug 30 10:41 compareSnaps
-rwxr-xr-x   1 fnrizzi  staff   348K Aug 30 10:42 computeThinSVD
-rwxr-xr-x   1 fnrizzi  staff   597K Aug 30 10:41 extractStateFromSnaps
-rwxr-xr-x   1 fnrizzi  staff   866K Aug 30 10:41 reconstructFomState
-rwxr-xr-x   1 fnrizzi  staff   627K Aug 30 10:41 reconstructSeismogram
-rwxr-xr-x   1 fnrizzi  staff   2.0M Aug 30 10:42 shwave_fom
-rwxr-xr-x   1 fnrizzi  staff   2.4M Aug 30 10:42 shwave_rom
drwxr-xr-x  11 fnrizzi  staff   352B Aug 30 10:41 tests
```
You can then run the tests to see if things are working as follows:
```bash
cd ${MYWORKDIR}/build
ctest
```



# Creating and running a full case
Here we describe how to create and run a case. There are a few steps involved,
so we will describe each one.

Define an env var as follows:
```bash
export MYRUNDIR=${MYWORKDIR}/myFirstRun
mkdir ${MYRUNDIR}
```

## Generating the mesh
The code has been developed such that the mesh is generated with Python
and is used by the C++ code. There are two main reasons for this choice:
first, it allows us to decouple the mesh generation from the actual physics code;
second, we developed the code to support the concept of sample mesh,
which is a key feature for nonlinear ROMs. This is **not** needed right now
to solve the current elastic shear wave problem because this is a linear problem,
but it can be useful if, in the fugute, we extend the code to support nonlinear problems.

To specify the grid, one only needs to specify the grid for the velocity points because
the stress points are defined based on the staggered scheme (see paper).
Assume the test case you want to run uses a grid of 150 x 600 velocity points
along the radial and polar directions, respectively.
To generate the mesh files proceed as follows:
```python
cd ${ESWSRCDIR}/meshing
python create_single_mesh.py -nr 150 -nth 600 -working-dir ${MYRUNDIR}
```
This should generate inside `$MYRUNDIR}` a directory called `mesh150x600`
containing the following files;
```bash
-rw-r--r--  1 fnrizzi  staff   2.0M Aug 30 11:19 coeff_vp.dat
-rw-r--r--  1 fnrizzi  staff    12M Aug 30 11:19 graph_sp.dat
-rw-r--r--  1 fnrizzi  staff   7.2M Aug 30 11:19 graph_vp.dat
-rw-r--r--  1 fnrizzi  staff   229B Aug 30 11:19 mesh_info.dat
```


<!-- - Second, look at -->
