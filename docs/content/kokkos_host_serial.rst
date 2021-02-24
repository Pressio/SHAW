Serial Host-only Kokkod Build
#############################

:breadcrumb: {filename}/kokkos_host_serial.rst
:date: 2021-02-12 11:00
:summary: Build the code for host serial-only execution using serial Kokkos backend
:category: building

###################
`1. Prerequisites`_
###################

* CMake>=3.13.0

* C++ (with support for c++14) compiler: we have tested this with GCC 8.3.1 and GCC 8.4.0

* BLAS/LAPACK: if you don't have them, we provide a script to build them for you

#########################
`2. Prepare environment`_
#########################

.. code:: bash

   export CXX=<path-to-your-C++-compiler>
   export ESWSRCDIR=<path-to-where-you-cloned-the-repository>
   export MYWORKDIR=<path-to-where-you-want-to-work-in> #e.g. ${HOME}/myWaveTest
   mkdir -p ${MYWORKDIR}

#################
`3. BLAS/LAPACK`_
#################

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

   export CC=<path-to-your-C-compiler>
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

Now, do and move to section 4:

.. code:: bash

   export BLAS_ROOT=${MYWORKDIR}/tpls/openblas/install
   export LAPACK_ROOT=${MYWORKDIR}/tpls/openblas/install
   export BLASLIBNAME=openblas
   export LAPACKLIBNAME=openblas


##############################
`4. Build Kokkos and Kernels`_
##############################

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
  you need to modify the flags passed to Kokkos inside
  `build_kokkos_and_kernels.sh <https://github.com/fnrizzi/SHAW/tree/master/bash_scripts/build_kokkos_and_kernels.sh>`_
  and rerun it.


#############################################
`5. Build the Shear Wave Code and Run Tests`_
#############################################

.. code:: bash

   cd ${ESWSRCDIR}/bash_scripts
   ./do_build.sh --working-dir=${MYWORKDIR} --kokkos-pfx=${KOKKOSPFX} --kokkos-ker-pfx=${KOKKOSKERPFX}
   cd ${MYWORKDIR}/build
   ctest

which should display (at the time of this writing we have these tests):

.. code:: bash


	  Start  1: meshinfo
     1/25 Test  #1: meshinfo ............................   Passed    0.26 sec
	  Start  2: parser_test_1
     2/25 Test  #2: parser_test_1 .......................   Passed    0.17 sec
	  Start  3: parser_test_2
     3/25 Test  #3: parser_test_2 .......................   Passed    0.16 sec
	  Start  4: parser_test_3
     4/25 Test  #4: parser_test_3 .......................   Passed    0.16 sec
	  Start  5: parser_test_4
     5/25 Test  #5: parser_test_4 .......................   Passed    0.16 sec
	  Start  6: parser_test_5
     6/25 Test  #6: parser_test_5 .......................   Passed    0.16 sec
	  Start  7: parser_test_6
     7/25 Test  #7: parser_test_6 .......................   Passed    0.16 sec
	  Start  8: seismogram_test
     8/25 Test  #8: seismogram_test .....................   Passed    0.17 sec
	  Start  9: forcing_rank1
     9/25 Test  #9: forcing_rank1 .......................   Passed    0.17 sec
	  Start 10: graphs
    10/25 Test #10: graphs ..............................   Passed    0.17 sec
	  Start 11: coords
    11/25 Test #11: coords ..............................   Passed    0.17 sec
	  Start 12: jacobian_vp
    12/25 Test #12: jacobian_vp .........................   Passed    0.17 sec
	  Start 13: jacobian_sp
    13/25 Test #13: jacobian_sp .........................   Passed    0.16 sec
	  Start 14: stress_labels
    14/25 Test #14: stress_labels .......................   Passed    0.17 sec
	  Start 15: fomInnerDomain
    15/25 Test #15: fomInnerDomain ......................   Passed    1.40 sec
	  Start 16: fomNearSurface
    16/25 Test #16: fomNearSurface ......................   Passed    1.28 sec
	  Start 17: fomNearCmb
    17/25 Test #17: fomNearCmb ..........................   Passed    1.68 sec
	  Start 18: fomSymmetryAxisThetaZero
    18/25 Test #18: fomSymmetryAxisThetaZero ............   Passed    2.26 sec
	  Start 19: fomSymmetryAxisThetaPi
    19/25 Test #19: fomSymmetryAxisThetaPi ..............   Passed    2.31 sec
	  Start 20: multiDepthsForcingRank1
    20/25 Test #20: multiDepthsForcingRank1 .............   Passed    1.95 sec
	  Start 21: multiPeriodsForcingRank1
    21/25 Test #21: multiPeriodsForcingRank1 ............   Passed    1.78 sec
	  Start 22: multiDepthsAndPeriodsForcingRank1
    22/25 Test #22: multiDepthsAndPeriodsForcingRank1 ...   Passed    5.19 sec
	  Start 23: multiDepthsForcingRank2
    23/25 Test #23: multiDepthsForcingRank2 .............   Passed    1.86 sec
	  Start 24: multiPeriodsForcingRank2
    24/25 Test #24: multiPeriodsForcingRank2 ............   Passed    1.03 sec
	  Start 25: multiDepthsAndPeriodsForcingRank2
    25/25 Test #25: multiDepthsAndPeriodsForcingRank2 ...   Passed    2.94 sec

   100% tests passed, 0 tests failed out of 25
