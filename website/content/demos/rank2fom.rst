Multi-Forcing Run
#################

:breadcrumb: {filename}/rank2fom.rst
:summary: Multi-Forcing FOM

.. container::

   This step-by-step demos shows how to run a simulation involving multiple forcings simultaneously.

   Here we assume you already followed the `step-by-step guide <{filename}/getstarted/build_kokkos_host_serial.rst>`_
   to build the code and used ``MYWORKDIR`` as the working directory for that procedure,
   so that ``${MYWORKDIR}/build`` contains all the executables.

`1. Prepare environment`_
=========================

.. code:: bash

   # ensure this env var points to the build directory
   export MYWORKDIR=<the-same-work-directory-used-for-building-process>

   export ESWSRCDIR=<path-to-the-source-code-repository>
   export MYRUNDIR=${MYWORKDIR}/myFirstRun
   mkdir -p ${MYRUNDIR}


`2. Generating the mesh`_
=========================

The code has been developed such that the mesh is generated with Python
and is used by the C++ code. There are two main reasons for this choice:
first, it allows us to decouple the mesh generation from the actual physics code;
second, we developed the code with an eye to later on study how to apply sample mesh
to nonlinear wav problems. This is **not** needed right now
to solve the current elastic shear wave problem because this is a linear problem,
but it can be useful if, in the future, we extend the code to support nonlinear problems.

To specify the grid, one only needs to specify the grid for the velocity points because
the stress points are defined based on the staggered scheme (see paper).
For this demo, we use a grid of ``200`` x ``1000`` velocity points
along the radial and polar directions, respectively.
To generate the mesh files proceed as follows:

.. code:: bash

   cd ${ESWSRCDIR}/meshing
   python create_single_mesh.py -nr 200 -nth 1000 -working-dir ${MYRUNDIR}


This should generate a directory ``${MYRUNDIR}/mesh200x1000`` containing:

.. code:: bash

   -rw-r--r--  1 fnrizzi  staff   4.5M Aug 30 12:20 coeff_vp.dat
   -rw-r--r--  1 fnrizzi  staff    28M Aug 30 12:20 graph_sp.dat
   -rw-r--r--  1 fnrizzi  staff    16M Aug 30 12:20 graph_vp.dat
   -rw-r--r--  1 fnrizzi  staff   231B Aug 30 12:20 mesh_info.dat


`3. Input file`_
================

todo
