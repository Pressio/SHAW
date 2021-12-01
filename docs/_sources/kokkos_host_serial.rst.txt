Host-only Serial
================

Prerequisites
--------------

* The SHAW repo: https://github.com/Pressio/SHAW

* CMake>=3.13.0

* C++ (with support for c++14) compiler: we have tested this with GCC 8.3.1 and GCC 8.4.0

* BLAS/LAPACK: if you don't have them, we provide a script to build them for you


Prepare environment
-------------------

.. code-block:: shell

   export CXX=<path-to-your-C++-compiler>
   export ESWSRCDIR=<path-to-where-you-cloned-the-repository>
   export MYWORKDIR=<path-to-where-you-want-to-work-in> #e.g. ${HOME}/myWaveTest
   mkdir -p ${MYWORKDIR}

Build the Shear Wave Code and Run Tests
---------------------------------------

.. code-block:: shell

   cd ${ESWSRCDIR}/bash_scripts
   ./do_build.sh --working-dir=${MYWORKDIR} --kokkos-pfx=${KOKKOSPFX} --kokkos-ker-pfx=${KOKKOSKERPFX}
   cd ${MYWORKDIR}/build
   ctest

which should display (at the time of this writing we have these tests):

.. code-block:: shell

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
