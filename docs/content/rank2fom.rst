Multi-forcing Run with rank-2
#############################

:breadcrumb: {filename}/rank2fom.rst
:summary: This demo shows how to simulate multiple trajectories using the rank-2 formulation.
:date: 2021-02-24 20:00
:category: demos

.. container::

   Here we assume you completed the building process, e.g., `step-by-step guide <{filename}/kokkos_host_serial.rst>`_,
   ``MYWORKDIR`` points to the working directory for that procedure, so that ``${MYWORKDIR}/build`` contains all the executables.

|

.. block-info:: Scope

		For the sake of demonstration, in this demo we solve *the same problem*
		solved in `the rank-1 demo <{filename}/rank1fommulti.rst>`_, except that
		we use the rank-2 formulation, which allows us to simulate all four sample
		trajectories in one single run.

`1. Prepare environment`_
=========================

.. code:: bash

   # ensure this env var points to the build directory
   export MYWORKDIR=<the-same-work-directory-used-for-building-process>
   export ESWSRCDIR=<path-to-the-source-code-repository>
   export MYRUNDIR=${MYWORKDIR}/mySecondRun
   mkdir -p ${MYRUNDIR}


`2. Generating the mesh`_
=========================
This is identical to `the mesh used in this demo <{filename}/rank1fommulti.rst>`_.


`3. Input file`_
================
The input file is identical to `the input for the rank-1 demo <{filename}/rank1fommulti.rst>`_,
except for one addition to the source yaml node:

.. code:: yaml

  #
  # general, io, material: as in the other demo
  #
  source:
    signal:
      # kind, depth, period, delay: same as the other one
      # ...

      # forcingSize defines many simultaneous trajectories to compute
      forcingSize: 4

The full input file can be copied:

.. code:: bash

   cp ${ESWSRCDIR}/demos/fom_rank2_sample_depth/input.yaml ${MYRUNDIR}


`3. Run the simulation`_
========================

.. code:: bash

   cd ${MYRUNDIR}
   ln -s ${MYWORKDIR}/build/shawExe .

   # if you use OpenMP build, remember to set
   # OMP_NUM_THREADS=4 OMP_PLACES=threads OMP_PROC_BIND=spread
   ./shawExe input.yaml
