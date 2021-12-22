Building: "expert" mode
=======================

Prerequisites
-------------

* This repo: ``git clone https://github.com/Pressio/SHAW``

* C++14 compiler: we have tested this with GCC 8.3.1, GCC 8.4.0, GCC 10.2.0.

* ``CMake>=3.16.0``

* ``BLAS`` and ``LAPACK``

* `Kokkos <https://github.com/kokkos/kokkos>`_ and
  `Kokkos-kernels <https://github.com/kokkos/kokkos-kernels>`_: last tested version ``3.5.00``

* `yaml-cpp <https://github.com/jbeder/yaml-cpp>`_: last tested version ``0.7.0``


Build
-----

.. code-block:: shell

   cmake \
   -DCMAKE_CXX_COMPILER=<fullpath-to-your-C++-compiler> \
   -DKOKKOSKERNELS_DIR=<fullpath-to-your-kernels-install-path> \
   -DYAMLCPP_DIR=<fullpath-to-your-yamlcpp-install-path> \
   -B <fullpath-to-where-you-want-to-build-the-code> \
   -S <fullpath-to-your-shaw-repository>

   # from within your build dir
   make -j4

   # running the tests is advised
   ctest
