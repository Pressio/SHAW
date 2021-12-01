Building: step-by-step
======================

If you are reading this page, it likely is because you want
a simplified (automated) way the get this done,
so that you can minimize the extra effort
in building the TPLs needed and the SHAW code.
This page tries to address this: it provides a step-by-step
guide that leverages some scripts we have prepared to simplify this.


Prerequisites
-------------

* This repo: ``git clone https://github.com/Pressio/SHAW``

* C++14 compiler with OpenMP support.

  - We have tested this with GCC 8.3.1, GCC 8.4.0, GCC 10.2.0.

* ``CMake>=3.16.0``: easily installed via package manager

* ``BLAS`` and ``LAPACK``

  * should already be present on your machine, if not, use your package manager


Step 1: Prepare environment
----------------------------

Let's make things easy:

.. code-block:: shell

   export CXX=<path-to-your-C++14-compiler>
   export SHAWDIR=<path-to-where-you-cloned-the-SHAW-repository>

   export WORKDIR=${HOME}/myFirstShawBuild
   mkdir -p ${WORKDIR}


Step 2: Build TPLs
--------------------------------

To simplify this part, we have prepared script that
automates getting the TPLs:

.. code-block:: shell

   cd ${SHAWDIR}/bash_scripts
   bash build_tpls.sh ${WORKDIR} openmp

This script will fetch, build and install inside ``WORKDIR/tpls``
all TPLs needed: Kokkos-core, Kokkos-kernelas and yaml-cpp.

.. Attention::

   This will build Kokkos for host-only use with the OpenMP backend
   but **without** any architecture specifications. This is on purpose,
   because this step is meant to be as generic and simple as possible to get
   you started quickly. If you want to customize things, read
   more on the `Kokkos github <https://github.com/kokkos>`_.


Once the build terminates, you should see
the following structure inside your ``WORKDIR``:

.. code-block::

  tree -c -L 2
  .
  └── tpls
      ├── kokkos-3.5.00
      ├── kokkos-3.5.00.tar.gz
      ├── kokkos-kernels-3.5.00
      ├── kokkos-kernels-3.5.00.tar.gz
      ├── kokkos-build
      ├── kokkos-install
      ├── kokkos-kernels-build
      ├── kokkos-kernels-install
      ├── yaml-cpp-0.7.0.tar.gz
      ├── yaml-cpp-yaml-cpp-0.7.0
      ├── yamlcpp-build
      └── yamlcpp-install


Step 3: Build SHAW
------------------

.. code-block:: shell

   cd ${WORKDIR}
   mkdir shaw-build && cd shaw-build

   cmake \
     -DKokkosKernels_DIR=${WORKDIR}/tpls/kokkos-kernels-install/lib/cmake/KokkosKernels/ \
     -Dyaml-cpp_DIR=${WORKDIR}/tpls/yamlcpp-install/share/cmake/ \
     ${SHAWDIR}

   make -j4

   # running the SHAW tests is advised
   ctest
