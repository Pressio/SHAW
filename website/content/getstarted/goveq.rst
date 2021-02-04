Governing Equations and Discretization
######################################

:breadcrumb: {filename}/goveq.rst
:summary: Equations and Discretization

.. role:: math-info(math)
    :class: m-default

.. container::

	*The objective is to model the axisymmetric evolution of elastic seismic shear waves.*

	First, we set the notation.
	Assuming the Earth (or another target planet) can be approximated as a sphere, the following
	figure shows our notation of the spherical coordinate system.

	.. figure:: {static}/img/sc.svg
		:width: 350 px


	In the axisymmetric approximation, one assumes fields/quantities to not vary along :math-info:`\phi`.
	This means that all the derivatives with respect to :math-info:`\phi` can be dropped.

	With this assumption, the set of equations governing the time evolution
	of elastic waves in the velocity-stress formulation can be written as:

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


	where :math-info:`t` represents time, :math-info:`r \in [0, r_{earth}]` is the radial distance bounded by the radius of the earth,
	:math-info:`\theta \in [0, \pi]` is the polar angle, :math-info:`\rho(r, \theta)` is the density,
	:math-info:`v(r, \theta, t)` is the velocity, :math-info:`\sigma_{r\phi}(r, \theta, t)` and
	:math-info:`\sigma_{\theta\phi}(r, \theta, t)` are the two components of the stress tensor remaining after the
	axisymmetric approximation, :math-info:`f(r, \theta,t)` is the forcing term,
	and the shear modulus is :math-info:`G(r, \theta) = v_s^2(r, \theta) \rho(r, \theta)`,
	with :math-info:`v_s` being the shear wave velocity.

	In practice, the axisymmetric approximation means that one solves the
	above governing equations over a *circular sector*.
	Such a formulation is referred to as 2.5-dimensional because it involves
	a 2-dimensional spatial domain (a circular sector of the Earth)
	but models point sources with correct 3-dimensional spreading {cite}.

	.. Note that we assume both the density and shear modulus to only depend on the spatial coordinates.

	Shear waves cannot propagate in liquids.
	Therefore, when modeling the Earth, the system of equations above is not
	applicable to the core region of the Earth, and are thus solved
	in the region bounded between the core-mantle boundary (CMB) located
	at :math-info:`r_{cmb} = 3,480` km,
	and the Earth surface located at :math-info:`r_{earth} = 6,371` km.

        The 2D spatial domain is discretized using a staggered grid as shown in the figure below.
	Staggered grids are typical for seismic modeling and wave problems in general.
	We use a second-order centered finite-difference method for the spatial operators.

	.. figure:: {static}/img/mesh.png
		:width: 450 px
