.. raw:: html

    <style> .red {color:red} </style>

.. role:: red


SHeAr Waves (SHAW) Simulator
============================

This project contributes an open-source C++ code using the
Kokkos programming model to simulate the generation and propagation
of elastic shear waves in an axi-symmetric domain.

.. figure:: ../img/logo1.png
   :align: center
   :width: 26%

   Sample countour plot of a velocity field obtained using the SHAW code.


Motivation
----------

Seismic modeling and simulation is an active field of research
because of its critical importance to understand the generation,
propagation and effects of seismic events (aka earthquakes on Earth,
moonquakes on the moon, etc), and artificial explosions.

Broadly, one can distinguish between two main types
of seismic waves: shear and pressure.
Shear waves are also called S-waves (or secondary) because they come
after P-waves (or primary). The main difference is that S-waves
are *transversal* (particles oscillate perpendicularly to the direction
of wave propagation), while P-waves are *longitudinal* (particles oscillate
in the same direction as the wave). Both P- and S-waves
are body waves, because they travel through the interior of the Earth
(or some other planet), and their evolution is affected
by the generating source as well as the material properties of the medium,
namely density, stiffness, composition, etc.

Modeling and simulating these systems is challenging because (a) physical models
contain a large number of parameters (e.g., anisotropic material properties,
signal forms and parametrizations); and (b) simulating at global scale
with high-accuracy requires a large computational cost.

We hope our code can help advance this field and foster new research and related work.


Highlights and features
-----------------------

*  Implementation based on the `Kokkos programming model <https://github.com/kokkos>`_
   for performance portability

*  :doc:`Velocity-stress formulation in an axi-symmetric domain <goveq>`

*  Support for the following material models:

   - :ref:`single layer model <singlelayerdescription>`

   - :ref:`bilayer model <twolayerdescription>`

   - :ref:`the Preliminary Reference Earth Model (PREM) <premdescription>`

   - :ref:`user-defined/custom model <customdescription>`

   These are 1D models because they only depend on the radial distance.
   The modularity of the code allows one to easily add new models

*  Simulating the dynamics in another planet/axisymmetric body is relatively easy:
   you have to create a mesh suitable for that planet, and a suitable material model

*  The code implements what we refer to as "rank-1" and "rank-2" formulations:

   *  *rank-1*:

      * :ref:`discrete state and forcing are stored as 1D arrays <rank1fom>`

      * this is useful to simulate the wave dynamics due to a *single forcing term*

      * :doc:`See the demo! <demo1>`

   *  *rank-2*:

      * :ref:`discrete state and forcing are stored using rank-2 tensors (i.e. matrices) <rank2fom>`

      * this is useful to *simultaneously* solve the wave
	dynamics for *multiple forcing realizations* (e.g. multiple
	source locations and/or periods). This rank-2 formulation
	has an advantage from a computational standpoint because
	it has higher computational intensity, thus benefiting
	efficient ensemble propagation

      * :doc:`See the demo! <demo3>`

How to cite
-----------

If you use this code, please cite the github page and the following paper:

.. code-block:: bibtex

   @article{RIZZI2021113973,
     title = {A compute-bound formulation of Galerkin model reduction for linear time-invariant dynamical systems},
     journal = {Computer Methods in Applied Mechanics and Engineering},
     volume = {384},
     pages = {113973},
     year = {2021},
     issn = {0045-7825},
     doi = {https://doi.org/10.1016/j.cma.2021.113973},
     url = {https://www.sciencedirect.com/science/article/pii/S0045782521003042},
     author = {Francesco Rizzi and Eric J. Parish and Patrick J. Blonigan and John Tencer}
   }


Contents
========

.. toctree::
    :maxdepth: 2

    goveq
    build_expert
    build_stepbystep
    inputfile
    demos
    performance
    GitHub Repo <https://github.com/Pressio/SHAW>
    Open an issue/feature req. <https://github.com/Pressio/SHAW/issues>
    license
