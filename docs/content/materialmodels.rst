Material Models
###############

:breadcrumb: {filename}/materialmodels.rst
:summary: Material Models: overview and supported parametrizations
:date: 2021-02-22 11:00

.. role:: math-info(math)
    :class: m-default

The code currently supports the following material models:

* single layer

* two layers

* PREM

* custom material model

.. block-danger:: No default choice is set.

		  When you create an input file for a simulation, you *need to* define a material model.



##################################
`1. Single Layer Material Model`_
##################################

A single material *without* discontinuities as shown in the
figure below, and such that both the density and shear velocity only have radial dependence.


.. figure:: {static}/img/mat_f1.png
	    :scale: 40 %


The code currently supports up to quadratic parametrizations as:

.. math::

   \rho(x) = a_0 + a_1 x + a_2 x^2

.. math::

   v_s(x) = b_0 + b_1 x + b_2 x^2


where the coefficients should be provided such that the density must be in *[kg/m^3]*,
and the shear velocity must be in *[m/s]*, considering :math-info:`x` to be in units of *[Km]*.

Defining such a model in the input file can be done as follows:

.. code:: yaml

  material:
    kind: unilayer
    layer:
      density: [a0, a1, a2]    # density  must have units of kg/m^3
      velocity: [b0, b1, b2]   # velocity must have units of m/s


Note that if you want a homogenouos material with constant density
and shear velocity, you can just set :math-info:`a_0, b_0` and set all other
coefficients to zero.

|

############################
`2. Bilayer Material Model`_
############################
Represents a material model with two layers, separated by a
single discontinuity as shown in the figure below.
Both the density and shear velocity only have radial dependence.


.. figure:: {static}/img/mat_f2.png
      :scale: 40 %


As above, the profiles can be up to quadratic, and the coefficients
should be provided such that the density is in *[kg/m^3]*,
the shear velocity in *[m/s]*, assuming :math-info:`x` to be in units of *[Km]*, and the
discontinuity is :math-info:`d` [km] deep.

Defining such a model in the input file can be done as follows:

.. code:: yaml

  material:
    kind: unilayer
    layer1:
      density: [a0, a1, a2]    # density  must have units of kg/m^3
      velocity: [b0, b1, b2]   # velocity must have units of m/s
    layer2:
      depth: d                 # must have units of km
      density: [c0, c1, c2]    # density  must have units of kg/m^3
      velocity: [d0, d1, d2]   # velocity must have units of m/s

where it is intended that within each layer, the density and shear velocity can
have different parametrizations. Note that this supports different
parametrizations within each layer, and potentially discontinuous profiles.


|

##############################################
`3. Preliminary Reference Earth Model (PREM)`_
##############################################
The PREM model is a radial model representing the average Earth properties, and one of the most
commonly adopted. Choose it from the input file as:

.. code:: yaml

  material:
    kind: prem

The details of the parametrization for the PREM model are `handled directly within the code <https://github.com/fnrizzi/SHAW/blob/master/src/shared/material_models/material_model_prem.hpp>`_.


For more details, check the following references:

* Dziewonski, A.M., and D.L. Anderson. 1981. “Preliminary reference Earth model.” Phys. Earth Plan. Int. 25:297-356.

* http://ds.iris.edu/ds/products/emc-prem/

* https://www.cfa.harvard.edu/~lzeng/papers/PREM.pdf

.. block-warning:: PREM is for Earth only

		   The PREM model only makes sense when you are simulating the Earth.
		   So your domain must be bounded between the core-mantle boundary (CMB)
		   located at :math-info:`r_{cmb} = 3,480` km and the Earth surface located at :math-info:`r_{earth} = 6,371` km.
		   These are the default bounds used by `the meshing script <https://github.com/fnrizzi/SHAW/blob/master/meshing/create_single_mesh.py>`_.


#############################################
`4. Using a Custom Model from the main file`_
#############################################

If you are *not* interested in using one of the choices described above, you can easily try a custom one without needing to change the source code. To do so, you need two do two things:

1.  in your input file, you need to set:

    .. code:: yaml

        material:
	  kind: custom


2.  and also modify the ``MyCustomMaterialModel`` inside `the main file <https://github.com/fnrizzi/SHAW/blob/master/src/kokkos/main_fom.cc>`_
    as you desire such that when the ``computeAt`` method is called for a given grid point in the domain, you set the local density and shear velocity according to you model.

.. block-info:: Extending the set of supported models

		The modular structure of the code allows to easily add new models: this can easily be done by adding a new
		derived class inside `the models <https://github.com/fnrizzi/SHAW/tree/master/src/shared/material_models>`_,
		add an ``enum`` field that identifies that model in `this file <https://github.com/fnrizzi/SHAW/blob/master/src/shared/enums/supported_material_model_enums.hpp>`_, and adding the code in `the parser class <https://github.com/fnrizzi/SHAW/blob/master/src/shared/parser/parser_material_model.hpp>`_ to recognize that if selected from the input file.
