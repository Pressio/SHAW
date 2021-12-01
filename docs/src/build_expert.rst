Building: "expert" mode
=======================

You need:

* This repo: ``git clone https://github.com/Pressio/SHAW``

* C++14 compiler: we have tested this with GCC 8.3.1 and GCC 8.4.0

* ``CMake>=3.16.0``

* ``BLAS`` and ``LAPACK``

* `Kokkos <https://github.com/kokkos/kokkos>`_ and
  `Kokkos-kernels <https://github.com/kokkos/kokkos-kernels>`_: last tested version ``3.5.00``

* `yaml-cpp <https://github.com/jbeder/yaml-cpp>`_: last tested version ``0.7.0``


Then, as usual for a typical CMake project, you can do:

.. code-block:: shell

   export CXX=<path-to-your-C++-compiler>
   export SHAWDIR=<path-to-where-you-cloned-the-shaw-repository>

   mkdir build && cd build

   cmake \
   -DKokkosKernels_DIR=<your-kernels-install-path>/lib/cmake/KokkosKernels/ \
   -Dyaml-cpp_DIR=<your-yamlcpp-install-path>/share/cmake/ \
   ${SHAWDIR}

   make -j4

   # running the tests is advised
   ctest


..
   export WORKDIR=<path-to-where-you-want-to-work-in>  #e.g. ${HOME}/myWaveTest
   mkdir -p ${WORKDIR}

   cd ${SHAWDIR}/bash_scripts
   ./do_build.sh --working-dir=${WORKDIR} --kokkos-pfx=${KOKKOSPFX} --kokkos-ker-pfx=${KOKKOSKERPFX}
   cd ${WORKDIR}/build
   ctest
