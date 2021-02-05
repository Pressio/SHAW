Elastic Shear Waves
###################

:save_as: index.html
:cover: {static}/img/top5.jpg
:url:
:description: Performance-portable simulator for elastic shear waves.
:summary: Performance-portable simulator for elastic shear waves.
:hide_navbar_brand: True
:landing:
    .. container:: m-row

        .. container:: m-col-l-9 m-push-l-1 m-nopadb

            .. raw:: html

                <h1 style="text-transform:capitalize">Elastic Shear Waves (ShWav)</h1>

    .. container:: m-row

        .. container:: m-col-l-8 m-push-l-1

            Seismic modeling and simulation is an active field of research
            because of its critical importance to understand the generation,
            propagation and effects of earthquakes and artificial explosions.

            One can distinguish between two main types of seismic waves: shear and pressure.
            Shear waves are also called S-waves (or secondary) because they come
            after P-waves (or primary). The main difference between them is that S-waves
            are *transversal* (particles oscillate perpendicularly to the direction
            of wave propagation), while P-waves are *longitudinal* (particles oscillate
            in the same direction as the wave). Both P- and S-waves
            are body waves, because they travel through the interior of the earth
            (or some other planet), and their evolution is affected
            by the generating source as well as the material properties of the medium,
            namely density, stiffness, composition, etc.

            Modeling and simulating these systems is challenging because (a) physical models
	    contain a large number of parameters (e.g., anisotropic material properties,
            signal forms and parametrizations); and (b) simulating these systems at global scale
            with high-accuracy requires a large computational cost.

            *This project contributes to the field by providing an open-source
            C++ code to simulate elastic shear waves in an axi-symmetric domain.*


        .. container:: m-col-l-3 m-push-l-1

            .. figure:: {static}/img/logo1.png
                        :scale: 50 %

    .. .. container:: m-row

    ..     .. container:: m-col-l-9 m-push-l-1

    ..         .. raw:: html

    ..             <p class="m-text m-default m-big"><i>This project presents an
    ..             open-source C++ code to simulate elastic shear waves in an axi-symmetric domain.</i></p>


    .. container:: m-row

        .. container:: m-col-l-11 m-push-l-1

            **HIGHLIGHTS AND CAPABILITIES:**

            * | The current implementation relies on Kokkos, but we are porting the code to other programming models.

            * The code implements what we call "rank-1" and "rank-2" formulations:

                * rank-1:
                    * the discrete state and forcing term are stored in 1D arrays
                    * this is used to simulate the case of a *single forcing realization*

                * rank-2:
                    * the discrete state and forcing term are stored in rank-2 tensors (i.e. matrices)
                    * this is useful to simultaneously solve the wave dyanmics
		      for *multiple forcing realizations* (e.g. multiple source locations and/or periods).
		      This rank-2 formulation has an advantage from a computational
                      standpoint because it has higher computational intensity,
		      thus benefiting efficient ensemble propagation of uncertainties

            * We use the velocity-shear formulation in an axi-symmetric domain, leading to a so-called 2.5 dimensional model.

            * The code currently supports the following radial material models: a single layer,
              a bilayer model and the Preliminary Reference Earth Model (PREM).
              We remark that these are 1D models in the sense that they only depend on the radial distance.
              Given the modularity of the code, one can easily add other material models.

	    * The modular strucuture of the code makes it easy to simulate
	      the wave dynamics in another planet/axisymmetric body:
	      one just has to create a grid suitable for that planet, and a suitable material model.

    .. container:: m-row

        .. container:: m-col-l-9 m-push-l-1

	    **GET STARTED:**

            1. read about the `governing equations, domain and discretization <{filename}/getstarted/goveq.rst>`_

            2. build the code using of the following versions:

		* `Host serial Kokkos <{filename}/getstarted/build_kokkos_host_serial.rst>`_:
		  Kokkos-only version with *host serial* backend.

		* `Host OpenMP Kokkos <{filename}/getstarted/build_kokkos_host_omp.rst>`_:
		  Kokkos-only version with *host OpenMP* backend.

            3. step-by-step demo for `single forcing (or rank-1) simulation <{filename}/demos/rank1fom.rst>`_

	    4. step-by-step demo for `multi-forcing (or rank-2) simulation <{filename}/demos/rank2fom.rst>`_

            5. learn about the structure of the input file


    .. container:: m-row

        .. container:: m-col-l-10 m-push-l-1

            If you use the code, please cite:
            *A compute-bound formulation of Galerkin model reduction for linear time-invariant dynamical systems*, by F.Rizzi, E.J.Parish, P.J.Blonigan, J.Tencer (https://arxiv.org/abs/2009.11742).
