Host Serial Kokkos Build
########################

:breadcrumb: {filename}/build_kokkos_host_serial.rst
:summary: Host Serial Kokkos

.. role:: math-info(math)
    :class: m-default


`1. Prerequisites`_
====================

* CMake>=3.13.0

* C, C++ (with support for c++14) compilers: we have tested this with GCC 8.3.1 and GCC 8.4.0

* BLAS/LAPACK: if you don't have them, we provide a script to build them for you

`2. Prep environment`_
======================

.. code:: bash

   export CC=<path-to-your-C-compiler>
   export CXX=<path-to-your-C++-compiler>

   export ESWSRCDIR=<path-to-where-you-clone-the-repository>

   export MYWORKDIR=<path-to-where-you-want-to-work-in> #e.g. ${HOME}/myWaveTest
   mkdir -p ${MYWORKDIR}


`3. BLAS/LAPACK`_
=================

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


`4: Build Kokkos and Kernels`_
==============================
Now that you BLAS/LAPACK is ready, we build Kokkos core and kernels as follows:

.. code:: bash

   cd ${MYWORKDIR}
   [[ ! -d tpls ]] && mkdir tpls
   cd tpls
   cp ${ESWSRCDIR}/bash_scripts/kokkos_host_serial/build_kokkos_and_kernels.sh .
   export KOKKOSPFX=${MYWORKDIR}/tpls/kokkos/kokkos_install
   export KOKKOSKERPFX=${MYWORKDIR}/tpls/kokkos/kokkos_kernels_install
   bash build_kokkos_and_kernels.sh

**Remarks**:

* the script above does a simple *serial build* to get you started quickly on any system.

* If you want to enable arch-specific optimizations following
  the `Kokkos userguide <https://github.com/kokkos/kokkos>`_
  and `here <https://github.com/kokkos/kokkos-kernels/wiki/Building>`_,
  you need to modify the flags passed to
  `build_kokkos_and_kernels.sh <https://github.com/fnrizzi/ElasticShearWaves/tree/master/bash_scripts/kokkos_host_serial/build_kokkos_and_kernels.sh>`_
  and rerun it.


`5: Build the Elastic Shear Wave Code`_
=======================================

.. code:: bash

   cd ${ESWSRCDIR}/bash_scripts
   ./do_build.sh --working-dir=${MYWORKDIR} --kokkos-pfx=${KOKKOSPFX} --kokkos-ker-pfx=${KOKKOSKERPFX}

this should generate inside ``${MYWORKDIR}/build`` the following:
