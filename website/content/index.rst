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

        .. container:: m-col-l-8 m-push-l-1 m-nopadb

            .. raw:: html

                <h1 style="text-transform:capitalize">Elastic Shear Waves Proxy App</h1>

    .. container:: m-row

        .. container:: m-col-l-7 m-push-l-1

            Seismic modeling and simulation is an active field of research
            because of its importance in understanding the generation,
            propagation and effects of earthquakes as well as artificial explosions.

            One can distinguish between two main types of seismic waves: shear and pressure.
            Shear waves are also called S-waves (or secondary) because they come
            after P-waves (or primary). The main difference between them is that S waves
            are transversal (particles oscillate perpendicularly to the direction
            of wave propagation), while the P waves are longitudinal (particles oscillate
            in the same direction as the wave). Both P and S waves
            are body waves, because they travel through the interior of the earth
            (or some other planet) and their trajectories are affected
            by the generating source as well as the material properties of the medium,
            namely density, stiffness, composition, etc.

            Modeling and simulating these systems is challenging because:
            (a) physical models contain a large number of parameters (e.g., anisotropic material properties,
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

        .. container:: m-col-l-10 m-push-l-1

            Highlights:

            * The current implementation relies on Kokkos, but we are porting the code to other programming models.<br/>
              For example, the first porting will be to a distributed memory using, e.g., Tpetra from Trilinos.

            * The code implements what we call "rank-1" and "rank-2" formulations:

                * rank-1:
                    * the discretized state variables are stored in a 1D array
                    * this formulation is used to simulate a problem for a *single forcing* term

                * rank-2:
                    * the discretized state variables are stored in a rank-2 tensor (i.e. a matrix)
                    * this formulation is useful to simulate the evolution for *multiple forcing*
                      terms simultaneously. The rank-2 formulation has an advantage from a computational
                      standpoint because it has higher computational intensity.

            * We use the velocity-shear formulation in an axi-symmetric domain.

	    * The model is so-called 2.5 dimensional model.

            * The code currently supports the following radial material models: a single layer,
              a bilayer model and the Preliminary Reference Earth Model (PREM).
              We remark that these are 1D models in the sense that they only depend on the radial distance.
              Given the modularity of the code, one can easily add other material models.

	    * The code can be used to simulate the wave dynamics in another planet:
	      this involves creating the grid and adding a suitable material model for that planet.

    .. container:: m-row

        .. container:: m-col-l-9 m-push-l-1

            Where to go from here:

            1. read about the `governing equations, domain and discretization <{filename}/getstarted/goveq.rst>`_

            2. to build the code, you currently have the following options:

		* `Host serial Kokkos <{filename}/getstarted/build_kokkos_host_serial.rst>`_:
		  to build the Kokkos version with *host serial* backend.

		* `Host OpenMP Kokkos <{filename}/getstarted/build_kokkos_host_omp.rst>`_:
		  to build the Kokkos version with *host OpenMP* backend.

            3. learn how to generate the grid for a simulation

            4. step-by-step guide on running a simulation

            5. learn about the structure of the input file


    .. container:: m-row

        .. container:: m-col-l-10 m-push-l-1

            If you use the code, please cite:
            *A compute-bound formulation of Galerkin model reduction for linear time-invariant dynamical systems*, by F.Rizzi, E.J.Parish, P.J.Blonigan, J.Tencer (https://arxiv.org/abs/2009.11742).
