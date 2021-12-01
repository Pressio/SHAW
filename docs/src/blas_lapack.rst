BLAS/LAPACK
===========

To handle BLAS/LAPACK, we envision the following three scenarios.

You want to let CMake figure it out
-----------------------------------

Just unset (to be safe) the following env vars, and then move to section 4.

.. code-block:: shell

   unset BLAS_ROOT
   unset LAPACK_ROOT
   unset BLASLIBNAME
   unset LAPACKLIBNAME

You have and want to use a specific BLAS/LAPACK
-----------------------------------------------

If you want a specific BLAS/LAPACK, set the following and then move to section 4.

.. code-block:: shell

   export BLAS_ROOT=<path-to-blas-install-dir>
   export LAPACK_ROOT=<path-to-lapack-install-dir>
   export BLASLIBNAME=<blas-lib-name-without-extension>     #e.g. openblas
   export LAPACKLIBNAME=<lapack-lib-name-without-extension> #e.g. lapack


You want to install them
------------------------

You can install BLAS/LAPACK yourself or you can follow the steps below
to install OpenBLAS. We chose OpenBLAS for simplicity
since it contains both BLAS and LAPACK and it is fairly easy to build.

.. code-block:: shell

   export CC=<path-to-your-C-compiler>
   export FC=<path-to-your-Fortran-compiler>
   cd ${MYWORKDIR}
   mkdir tpls && cd tpls
   cp ${ESWSRCDIR}/bash/build_openblas.sh .
   bash build_openblas.sh


If this succeeds, inside ``${MYWORKDIR}/tpls/openblas/install/lib``
you see something as:

.. code-block:: shell

   drwxr-xr-x  3 fnrizzi  staff    96B Aug 30 09:40 cmake
   lrwxr-xr-x  1 fnrizzi  staff    34B Aug 30 09:40 libopenblas.0.dylib
   lrwxr-xr-x  1 fnrizzi  staff    30B Aug 30 09:40 libopenblas.a
   lrwxr-xr-x  1 fnrizzi  staff    34B Aug 30 09:40 libopenblas.dylib
   drwxr-xr-x  3 fnrizzi  staff    96B Aug 30 09:40 pkgconfig


Then do: 

.. code-block:: shell

   export BLAS_ROOT=${MYWORKDIR}/tpls/openblas/install
   export LAPACK_ROOT=${MYWORKDIR}/tpls/openblas/install
   export BLASLIBNAME=openblas
   export LAPACKLIBNAME=openblas

