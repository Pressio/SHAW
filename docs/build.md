
# Step-by-step building process
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
The following has been tested with GCC-8.3.1 and GCC-8.4.0.

## Step 2: BLAS/LAPACK
Here we build and install BLAS and LAPACK.
Proceed as follows:
```bash
cd ${MYWORKDIR}
mkdir tpls && cd tpls
cp ${ESWSRCDIR}/bash/build_openblas.sh .
bash build_openblas.sh
```
this should build and install OpenBLAS such that
inside `${MYWORKDIR}/tpls/openblas/install/lib` you should see something as:
```bash
total 78328
drwxr-xr-x  3 fnrizzi  staff    96B Aug 30 09:40 cmake
lrwxr-xr-x  1 fnrizzi  staff    34B Aug 30 09:40 libopenblas.0.dylib
lrwxr-xr-x  1 fnrizzi  staff    30B Aug 30 09:40 libopenblas.a
lrwxr-xr-x  1 fnrizzi  staff    34B Aug 30 09:40 libopenblas.dylib
drwxr-xr-x  3 fnrizzi  staff    96B Aug 30 09:40 pkgconfig
```
and then set:
```bash
export BLAS_ROOT=${MYWORKDIR}/tpls/openblas/install
export BLASLIBNAME=openblas
export LAPACK_ROOT=${MYWORKDIR}/tpls/openblas/install
export LAPACKLIBNAME=openblas
```
**Note**: here, we use openblas for simplicity but if you already
have BLAS/LAPACK installed, you can skip the build step above
and directly set the needed env vars to
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
* for  now, we only enable the OpenMP backend;
* if you already have Kokkos/Kokkos-kernels installed, you can skip the build step
above and directly set the needed env vars to point to your installation.


## Step 4: Build the Elastic Shear Wave executables
Proceed as follows:
```bash
cd $ESWSRCDIR}
./do_build.sh --working-dir=${MYWORKDIR} --kokkos-pfx=${KOKKOSPFX} --kokkos-ker-pfx=${KOKKOSKERPFX} --omp=yes
```
this should generate inside `${MYWORKDIR}/build` the following:
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
which should display (if tests pass) the following:
```bash
      Start  1: fomInnerDomainKokkos1
 1/10 Test  #1: fomInnerDomainKokkos1 .............   Passed    0.66 sec
      Start  2: fomInnerDomainKokkos2
 2/10 Test  #2: fomInnerDomainKokkos2 .............   Passed    0.60 sec
      Start  3: fomNearSurfaceKokkos1
 3/10 Test  #3: fomNearSurfaceKokkos1 .............   Passed    0.67 sec
      Start  4: fomNearSurfaceKokkos2
 4/10 Test  #4: fomNearSurfaceKokkos2 .............   Passed    0.63 sec
      Start  5: fomNearCmbKokkos1
 5/10 Test  #5: fomNearCmbKokkos1 .................   Passed    0.91 sec
      Start  6: fomNearCmbKokkos2
 6/10 Test  #6: fomNearCmbKokkos2 .................   Passed    0.80 sec
      Start  7: fomSymmetryAxisThetaZeroKokkos1
 7/10 Test  #7: fomSymmetryAxisThetaZeroKokkos1 ...   Passed    1.14 sec
      Start  8: fomSymmetryAxisThetaZeroKokkos2
 8/10 Test  #8: fomSymmetryAxisThetaZeroKokkos2 ...   Passed    1.09 sec
      Start  9: fomSymmetryAxisThetaPiKokkos1
 9/10 Test  #9: fomSymmetryAxisThetaPiKokkos1 .....   Passed    1.17 sec
      Start 10: fomSymmetryAxisThetaPiKokkos2
10/10 Test #10: fomSymmetryAxisThetaPiKokkos2 .....   Passed    1.08 sec
```
