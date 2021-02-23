Host Serial Kokkos Build
########################

:breadcrumb: {filename}/kokkos_host_serial.rst
:summary: Building with Host Serial Kokkos
:date: 2021-02-12 11:00


`1. Prerequisites`_
====================

* CMake>=3.13.0

* C, C++ (with support for c++14) compilers: we have tested this with GCC 8.3.1 and GCC 8.4.0

* BLAS/LAPACK: if you don't have them, we provide a script to build them for you


`2. Prepare environment`_
=========================

.. code:: bash

   export CC=<path-to-your-C-compiler>
   export CXX=<path-to-your-C++-compiler>

   export ESWSRCDIR=<path-to-where-you-cloned-the-repository>

   export MYWORKDIR=<path-to-where-you-want-to-work-in> #e.g. ${HOME}/myWaveTest
   mkdir -p ${MYWORKDIR}


`3. BLAS/LAPACK`_
=================

To handle BLAS/LAPACK, we envision the following three scenarios.

`a. You want to let CMake figure it out`_
-----------------------------------------
Just unset (to be safe) the following env vars, and then move to section 4.

.. code:: bash

   unset BLAS_ROOT
   unset LAPACK_ROOT
   unset BLASLIBNAME
   unset LAPACKLIBNAME

`b. You have and want to use a specific BLAS/LAPACK`_
-----------------------------------------------------
If you want a specific BLAS/LAPACK, set the following and then move to section 4.

.. code:: bash

   export BLAS_ROOT=<path-to-blas-install-dir>
   export LAPACK_ROOT=<path-to-lapack-install-dir>
   export BLASLIBNAME=<blas-lib-name-without-extension>     #e.g. openblas
   export LAPACKLIBNAME=<lapack-lib-name-without-extension> #e.g. lapack


`c. You want to install them`_
------------------------------

You can install BLAS/LAPACK yourself or you can follow the steps below
to install OpenBLAS. We chose OpenBLAS for simplicity
since it contains both BLAS and LAPACK and it is fairly easy to build.

.. code:: bash

   export FC=<path-to-your-Fortran-compiler>
   cd ${MYWORKDIR}
   mkdir tpls && cd tpls
   cp ${ESWSRCDIR}/bash/build_openblas.sh .
   bash build_openblas.sh


If this succeeds, inside ``${MYWORKDIR}/tpls/openblas/install/lib``
you see something as:

.. code:: bash

   drwxr-xr-x  3 fnrizzi  staff    96B Aug 30 09:40 cmake
   lrwxr-xr-x  1 fnrizzi  staff    34B Aug 30 09:40 libopenblas.0.dylib
   lrwxr-xr-x  1 fnrizzi  staff    30B Aug 30 09:40 libopenblas.a
   lrwxr-xr-x  1 fnrizzi  staff    34B Aug 30 09:40 libopenblas.dylib
   drwxr-xr-x  3 fnrizzi  staff    96B Aug 30 09:40 pkgconfig

And then do:

.. code:: bash

   export BLAS_ROOT=${MYWORKDIR}/tpls/openblas/install
   export LAPACK_ROOT=${MYWORKDIR}/tpls/openblas/install
   export BLASLIBNAME=openblas
   export LAPACKLIBNAME=openblas


`4. Build Kokkos and Kernels`_
==============================
Now that you BLAS/LAPACK is ready, we build Kokkos core and kernels as follows:

.. code:: bash

   cd ${MYWORKDIR}
   [[ ! -d tpls ]] && mkdir tpls
   cd tpls
   cp ${ESWSRCDIR}/bash_scripts/build_kokkos_and_kernels.sh .
   export KOKKOSPFX=${MYWORKDIR}/tpls/kokkos/kokkos_install
   export KOKKOSKERPFX=${MYWORKDIR}/tpls/kokkos/kokkos_kernels_install
   bash build_kokkos_and_kernels.sh serial

**Remarks**:

* the script above does a simple *serial build* to get you started quickly on any system.

* If you want to enable arch-specific optimizations following
  the `Kokkos userguide <https://github.com/kokkos/kokkos>`_
  and `here <https://github.com/kokkos/kokkos-kernels/wiki/Building>`_,
  you need to modify the flags passed to
  `build_kokkos_and_kernels.sh <https://github.com/fnrizzi/ElasticShearWaves/tree/master/bash_scripts/build_kokkos_and_kernels.sh>`_
  and rerun it.


`5. Build the Shear Wave Code and Run Tests`_
=============================================

.. code:: bash

   cd ${ESWSRCDIR}/bash_scripts
   ./do_build.sh --working-dir=${MYWORKDIR} --kokkos-pfx=${KOKKOSPFX} --kokkos-ker-pfx=${KOKKOSKERPFX}
   cd ${MYWORKDIR}/build
   ctest

which should display (at the time of this writing we have these tests):

.. code:: bash

   Start  1: parser_test_1
   1/21 Test  #1: parser_test_1 .....................   Passed    0.32 sec
   Start  2: parser_test_2
   2/21 Test  #2: parser_test_2 .....................   Passed    0.19 sec
   Start  3: parser_test_3
   3/21 Test  #3: parser_test_3 .....................   Passed    0.22 sec
   Start  4: parser_test_4
   4/21 Test  #4: parser_test_4 .....................   Passed    0.19 sec
   Start  5: seismogram_test
   5/21 Test  #5: seismogram_test ...................   Passed    0.20 sec
   Start  6: forcing_rank1
   6/21 Test  #6: forcing_rank1 .....................   Passed    0.20 sec
   Start  7: graphs
   7/21 Test  #7: graphs ............................   Passed    0.19 sec
   Start  8: coords
   8/21 Test  #8: coords ............................   Passed    0.20 sec
   Start  9: jacobian_vp
   9/21 Test  #9: jacobian_vp .......................   Passed    0.20 sec
   Start 10: jacobian_sp
   10/21 Test #10: jacobian_sp .......................   Passed    0.20 sec
   Start 11: stress_labels
   11/21 Test #11: stress_labels .....................   Passed    0.20 sec
   Start 12: fomInnerDomainKokkos1
   12/21 Test #12: fomInnerDomainKokkos1 .............   Passed    0.67 sec
   Start 13: fomInnerDomainKokkos2
   13/21 Test #13: fomInnerDomainKokkos2 .............   Passed    0.47 sec
   Start 14: fomNearSurfaceKokkos1
   14/21 Test #14: fomNearSurfaceKokkos1 .............   Passed    0.50 sec
   Start 15: fomNearSurfaceKokkos2
   15/21 Test #15: fomNearSurfaceKokkos2 .............   Passed    0.47 sec
   Start 16: fomNearCmbKokkos1
   16/21 Test #16: fomNearCmbKokkos1 .................   Passed    0.64 sec
   Start 17: fomNearCmbKokkos2
   17/21 Test #17: fomNearCmbKokkos2 .................   Passed    0.64 sec
   Start 18: fomSymmetryAxisThetaZeroKokkos1
   18/21 Test #18: fomSymmetryAxisThetaZeroKokkos1 ...   Passed    0.86 sec
   Start 19: fomSymmetryAxisThetaZeroKokkos2
   19/21 Test #19: fomSymmetryAxisThetaZeroKokkos2 ...   Passed    0.85 sec
   Start 20: fomSymmetryAxisThetaPiKokkos1
   20/21 Test #20: fomSymmetryAxisThetaPiKokkos1 .....   Passed    0.84 sec
   Start 21: fomSymmetryAxisThetaPiKokkos2
   21/21 Test #21: fomSymmetryAxisThetaPiKokkos2 .....   Passed    0.85 sec

   100% tests passed, 0 tests failed out of 21
