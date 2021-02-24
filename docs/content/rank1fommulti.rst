Multi-forcing Run with rank-1
#############################

:breadcrumb: {filename}/rank1fommulti.rst
:summary: This demo shows how to simulate multiple trajectories using the rank-1 formulation.
:date: 2021-02-12 11:00
:category: demos

.. container::

   Here we assume you completed the building process, e.g., `step-by-step guide <{filename}/kokkos_host_serial.rst>`_,
   ``MYWORKDIR`` points to the working directory for that procedure, so that ``${MYWORKDIR}/build`` contains all the executables.

|

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

To specify the grid, one only needs to specify the grid for the velocity points because
the stress points are defined based on the staggered scheme (see paper).
For this demo, we use a grid of ``256`` x ``1024`` velocity points
along the radial and polar directions, respectively.
To generate the mesh files proceed as follows:

.. code:: bash

   cd ${ESWSRCDIR}/meshing
   python create_single_mesh.py -nr 256 -nth 1024 -working-dir ${MYRUNDIR}


This should generate a directory ``${MYRUNDIR}/mesh256x1024`` containing:

.. code:: bash

   .
   ├── [5.9M]  coeff_vp.dat
   ├── [ 37M]  graph_sp.dat
   ├── [ 21M]  graph_vp.dat
   └── [ 231]  mesh_info.dat


`3. Input file`_
================
For this demo, we use the following input file:

.. code:: yaml

  general:
    meshDir: ./mesh256x1024
    dt: 0.25
    finalTime: 2000.0
    checkNumericalDispersion: true
    checkCfl: true

  io:
    snapshotMatrix:
      binary: true
      velocity: {freq: 100, fileName: snaps_vp}
      stress:   {freq: 100, fileName: snaps_sp}

  seismogram:
    binary: false
    freq: 4
    receivers: [5,30,55,80,105,130,155,175]

  source:
    signal:
      kind: ricker

      # here we pass a list of depths to use as samples
      # this will automatically activate sampling
      depth: [240.,440.,540.,740.]

      period: 65.0
      delay: 180.0

  material:
    kind: prem

Which is ready to get:

.. code:: bash

   cp ${ESWSRCDIR}/demos/fom_rank1_sample_depth/input.yaml ${MYRUNDIR}


`3. Run the simulation`_
========================

.. code:: bash

   cd ${MYRUNDIR}
   ln -s ${MYWORKDIR}/build/shawExe .

   # if you use OpenMP build, remember to set
   # OMP_NUM_THREADS=4 OMP_PLACES=threads OMP_PROC_BIND=spread
   ./shawExe input.yaml


`5. Simulation data`_
=======================

After running the demo (have some patience because it takes some a couple minutes),
you should have inside ``${MYRUNDIR}`` the following files:

.. code:: bash

   coords_sp.txt #: coordinates of the velocity grid points
   coords_vp.txt #: coordinates of the stresses grid points

   seismogram_0  #: seismogram for depth = 240
   seismogram_1  #: seismogram for depth = 440
   seismogram_2  #: seismogram for depth = 540
   seismogram_3  #: seismogram for depth = 740

   snaps_vp_0    #: velocity snapshots for depth = 240
   snaps_vp_1    #: velocity snapshots for depth = 440
   snaps_vp_2    #: velocity snapshots for depth = 540
   snaps_vp_3    #: velocity snapshots for depth = 740

   snaps_sp_0    #: stresses snapshots for depth = 240
   snaps_sp_1    #: stresses snapshots for depth = 440
   snaps_sp_2    #: stresses snapshots for depth = 540
   snaps_sp_3    #: stresses snapshots for depth = 740


`4. Post-process data`_
=======================

To post-process the data, get the Python scripts created
for this demo and visualize the seismogram:

.. code:: bash

   cd ${MYRUNDIR}
   cp ${ESWSRCDIR}/demos/fom_rank1_sample_depth/plotSeismogram.py .
   python plotSeismogram.py


.. figure:: {static}/img/demo2_f1.png
