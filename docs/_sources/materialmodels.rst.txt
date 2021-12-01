Material Models
###############

The code currently supports the following material models:

* single layer

* two layers

* PREM

* custom material model


Single Layer Material Model
---------------------------

Models a domain with a single material *without*
discontinuities as shown in the figure below, and such that
both the density and shear velocity only have radial dependence.


.. image:: ../img/mat_f1.png
	   :scale: 40 %
	   :align: center


The code currently supports up to quadratic parametrizations as:

.. math::

   \rho(x) = a_0 + a_1 x + a_2 x^2 \\
   v_s(x) = b_0 + b_1 x + b_2 x^2

.. Important::
   the coefficients above should be provided considering the density
   must be in *[kg/m^3]*, and the shear velocity must be in *[m/s]*,
   and :math:`x` to be in units of *[Km]*.

Specifying such a model in the :doc:`input file <../inputfile>` simply involves specifing the coefficient as follows:

.. code-block:: yaml

  material:
    kind: unilayer
    layer:
      density:  [a0, a1, a2]   # density  must have units of kg/m^3
      velocity: [b0, b1, b2]   # velocity must have units of m/s

.. Tip::
   If you want a homogenouos material (i.e., constant density
   and shear velocity), you can just fix :math:`a_0, b_0`
   and set all other coefficients to zero.

|

Bilayer Material Model
----------------------

Represents a material model with two layers, separated by a
single discontinuity as shown in the figure below.
Both the density and shear velocity only have radial dependence.


.. image:: ../img/mat_f2.png
	   :scale: 40 %
	   :align: center


.. Important::
   Within each region, the profiles can be up to quadratic.
   The coefficients should be provided such that the density is in *[kg/m^3]*,
   the shear velocity in *[m/s]*, assuming :math:`x` to be in units of *[Km]*, and the
   discontinuity is :math:`d` [km] deep.

Defining such a model in the :doc:`input file <../inputfile>` can be done as follows:

.. code-block:: yaml

  material:
    kind: unilayer
    layer1:
      density:  [a0, a1, a2]   # density  must have units of kg/m^3
      velocity: [b0, b1, b2]   # velocity must have units of m/s
    layer2:
      depth:    integer        # must have units of km
      density:  [c0, c1, c2]   # density  must have units of kg/m^3
      velocity: [d0, d1, d2]   # velocity must have units of m/s

where it is intended that within each layer, the density and shear velocity can
have different parametrizations. Note that this supports different
parametrizations within each layer, and potentially discontinuous profiles.

|

Preliminary Reference Earth Model (PREM)
----------------------------------------

The PREM model is a radial model representing the average Earth properties, and one of the most
commonly adopted. Choose it from the input file as:

.. code-block:: yaml

  material:
    kind: prem

The details of the parametrization for the PREM model are `handled directly within the code <https://github.com/fnrizzi/SHAW/blob/master/src/shared/material_models/material_model_prem.hpp>`_.


For more details, check the following references:

* Dziewonski, A.M., and D.L. Anderson. 1981. “Preliminary reference Earth model.” Phys. Earth Plan. Int. 25:297-356.

* http://ds.iris.edu/ds/products/emc-prem/

* https://www.cfa.harvard.edu/~lzeng/papers/PREM.pdf


.. Caution::

   PREM is for Earth only!

   The PREM model only makes sense when you are simulating the Earth.
   So your domain must be bounded between the core-mantle boundary (CMB)
   located at :math:`r_{cmb} = 3,480` km and the Earth surface located at :math:`r_{earth} = 6,371` km.
   These are the default bounds used by `the meshing script <https://github.com/fnrizzi/SHAW/blob/master/meshing/create_single_mesh.py>`_.

|

Using a Custom Model
--------------------

If you are *not* interested in using one of the models above,
you can also try your own **without** needing to change the internal source code.
To do so, you need two do two things:

1.  in your input file, you need to set:

    .. code-block:: yaml

        material:
      	  kind: custom


2.  modify the ``MyCustomMaterialModel`` inside `the main file <https://github.com/fnrizzi/SHAW/blob/master/src/kokkos/main_fom.cc>`_
    as you desire such that when the ``computeAt`` method is called for a given grid point in the domain, you set the local density and shear velocity according to you model.


.. important::

  Extending the set of supported models

  The modular structure of the code allows to easily add new models: this can easily be done by adding a new
  derived class inside `the models <https://github.com/fnrizzi/SHAW/tree/master/src/shared/material_models>`_,
  add an ``enum`` field that identifies that model in `this file <https://github.com/fnrizzi/SHAW/blob/master/src/shared/enums/supported_material_model_enums.hpp>`_, and adding the code in `the parser class <https://github.com/fnrizzi/SHAW/blob/master/src/shared/parser/parser_material_model.hpp>`_ to recognize that if selected from the input file.
