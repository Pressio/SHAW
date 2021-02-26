OpenMP Host-only Kokkos Build
#############################

:breadcrumb: {filename}/kokkos_host_omp.rst
:date: 2021-02-12 11:00
:summary: Build the code for host-only execution using OpenMP Kokkos backend
:category: building

###################
`1. Prerequisites`_
###################

* The SHAW repo: https://github.com/Pressio/SHAW

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
This step is the same as described in `the serial host-only build <{filename}/kokkos_host_serial.rst>`_,

##############################
`4. Build Kokkos and Kernels`_
##############################

Now that you BLAS/LAPACK is ready, we build Kokkos core and kernels
with OpenMP support as follows:

.. code:: bash

   cd ${MYWORKDIR}
   [[ ! -d tpls ]] && mkdir tpls
   cd tpls
   cp ${ESWSRCDIR}/bash_scripts/build_kokkos_and_kernels.sh .
   export KOKKOSPFX=${MYWORKDIR}/tpls/kokkos/kokkos_install
   export KOKKOSKERPFX=${MYWORKDIR}/tpls/kokkos/kokkos_kernels_install
   bash build_kokkos_and_kernels.sh openmp

**Remarks**:

* If you want to enable arch-specific optimizations following
  the `Kokkos userguide <https://github.com/kokkos/kokkos>`_
  and `here <https://github.com/kokkos/kokkos-kernels/wiki/Building>`_,
  you need to modify the flags passed to Kokkos inside
  `build_kokkos_and_kernels.sh <https://github.com/Pressio/SHAW/tree/master/bash_scripts/build_kokkos_and_kernels.sh>`_
  and rerun it.

#############################################
`5. Build the Shear Wave Code and Run Tests`_
#############################################

This step is the same as described for `the serial host-only version <{filename}/kokkos_host_serial.rst>`_.
