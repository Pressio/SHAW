Governing Equations and Discretization
######################################

:breadcrumb: {filename}/goveq.rst
:summary: Equations and Discretization
:date: 2021-02-22 11:00

.. role:: yellow
    :class: m-text m-warning

.. role:: red
    :class: m-text m-danger

.. role:: blue
    :class: m-text m-info

.. role:: math-info(math)
    :class: m-default

.. container::

	:yellow:`Our objective is to model the evolution of elastic seismic shear waves in an axisymmetric domain.`

	Assuming the target body/planet (e.g. Earth) can be approximated as a sphere,
	we adopt a spherical coordinate system as shown in the figure below:

	.. figure:: {static}/img/sc.svg
		:width: 350 px


        In the axisymmetric approximation, one assumes that fields/quantities
	do not vary along :math-info:`\phi`, implying that all the derivatives
	with respect to :math-info:`\phi` can be dropped.

	With this assumption, the set of equations governing the time evolution
	of elastic waves in the velocity-stress formulation can be written as:

	.. block-info:: Governing Equations

			.. math::

				\rho (r, \theta) \frac{\partial v}{\partial t} (r, \theta,t) =
				\frac{\partial \sigma_{r\phi}}{\partial r}(r, \theta,t)
				+ \frac{\partial \sigma_{\theta\phi}}{r \partial \theta}(r, \theta,t)
				+ \frac{3}{r} \Big(\sigma_{r\phi}(r, \theta,t)
				+ 2 \sigma_{\theta\phi}(r, \theta,t) \cot{\theta} \Big) + f(r, \theta, t)

			.. math::

				\frac{\partial \sigma_{r\phi}}{\partial t}(r, \theta,t) =
				G(r, \theta)
				\left( \frac{\partial v}{\partial r}(r, \theta,t) - \frac{1}{r} v(r, \theta,t) \right)

			.. math::
				  \frac{\partial \sigma_{\theta\phi}}{\partial t}(r, \theta,t) =
				  G(r, \theta) \left( \frac{1}{r} \frac{\partial v}{\partial \theta}(r, \theta,t)
				  - \frac{\cot{\theta}}{r} v(r, \theta,t) \right)


			where:

 			   - :math-info:`t` represents time

			   - :math-info:`r \in [0, r_{surface}]` is the radial distance from origin to surface of the body

			   - :math-info:`\theta \in [0, \pi]` is the polar angle

			   - :math-info:`\rho(r, \theta)` is the density

			   - :math-info:`v(r, \theta, t)` is the velocity (for simplicity we drop the subscript,
			     but it is intended to be the :math-info:`v_{\phi}` velocity component)

			   - :math-info:`\sigma_{r\phi}(r, \theta, t)` and
			     :math-info:`\sigma_{\theta\phi}(r, \theta, t)` are the two components of the stress tensor remaining after the
			     axisymmetric approximation

			   - :math-info:`f(r, \theta,t)` is the forcing term

			   - :math-info:`G(r, \theta) = v_s^2(r, \theta) \rho(r, \theta)` is the shear modulus
			     and :math-info:`v_s` being the shear wave velocity.



	:yellow:`In practice, the axisymmetric approximation means that one solves the
	above governing equations over a *circular sector/block arc*.`
	Such a formulation is referred to as 2.5-dimensional because it involves
	a 2-dimensional spatial domain (a circular sector of the Earth)
	but models point sources with correct 3-dimensional spreading {cite}.

	.. Note that we assume both the density and shear modulus to only depend on the spatial coordinates.

	Shear waves cannot propagate in liquids.
	Therefore, when modeling the Earth, the system of equations above is not
	applicable to the core region of the Earth, and is solved in the region
	bounded between the core-mantle boundary (CMB) located at :math-info:`r_{cmb} = 3,480` km
	and the Earth surface located at :math-info:`r_{earth} = 6,371` km.
	When modeling a body that is not the Earth, one needs to use proper domain bounds.
	Note, however, that the code *always* assumes the domain NOT to include the
	singularity at the origin, which yields substantial simplifications.

        The 2D spatial domain is discretized using a staggered grid.
	Staggered grids are typical for seismic modeling and wave problems in general.
	We use a second-order centered finite-difference method for the spatial operators.
	Directly at the symmetry axis, i.e., :math-info:`\theta = 0, \pi`, the velocity
	is set to zero since it undefined here due to the cotangent term in its governing equation.
	This implies that at the symmetry axis the stress :math-info:`\sigma_{r,\phi}` is also zero.
	At the core-mantle boundary and earth surface we impose a free surface boundary
	condition (i.e., waves fully reflect), by setting the zero-stress condition
	:math-info:`\sigma_{r_{cmb},\phi} = \sigma_{r_{earth},\phi} = 0`.
	Note that this condition on the stress directly defines the velocity
	at the core-mantle boundary and earth surface and, therefore,
	no boundary condition on the velocity itself must be set there.
	We remark that, differently than (cite), we do not rely on ghost
	points to impose boundary conditions, but account for the boundary
	conditions directly when assembling the system matrix.
	As an example, the figure below shows the grid when modeling the Earth
	For the Earth, we know that there is a liquid core, so that region is considered.

	.. figure:: {static}/img/mesh.png
		    :width: 450 px

		    Schematic of the axi-symmetric domain for the Earth and staggered grid used for its discretization.
