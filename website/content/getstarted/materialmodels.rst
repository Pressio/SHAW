Material Models
###############

:breadcrumb: {filename}/materialmodels.rst
:summary: Material Models
:date: 2021-02-22 11:00

.. role:: math-info(math)
    :class: m-default

The code currently supports the following material models:

* single layer

* two layers

* PREM

* custom material model.


| When you create an input file to run a simulation, you need to define a material model.
| Note that the modular structure of the code allows for easily adding more models.


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
      density: [a0, a1, a2]    # must have units of kg/m^3
      velocity: [b0, b1, b2]   # must have units of m/s


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


As above, the profiles can be quadratic, and the coefficients should be provided
such that the density
must be in *[kg/m^3]*, and the shear velocity must be in *[m/s]*,
considering :math-info:`x` to be in units of *[Km]*, and the
discontinuity is located at :math-info:`d` [km] deep.

Define such a model in the input file can be done as follows:

.. code:: yaml

  material:
    kind: unilayer
    layer1:
      density: [a0, a1, a2]    # must have units of kg/m^3
      velocity: [b0, b1, b2]   # must have units of m/s
    layer2:
      depth: d           # must have units of km
      density: [c0, c1, c2]    # must have units of kg/m^3
      velocity: [d0, d1, d2]   # must have units of m/s

where it is intended that within each layer, the density and shear velocity can
have different parametrizations.


|

##############################################
`3. Preliminary Reference Earth Model (PREM)`_
##############################################
The PREM model is a radial model representing the average Earth properties, and one of the most
commonly adoptedo ones. For more details, check the following references:

* Dziewonski, A.M., and D.L. Anderson. 1981. “Preliminary reference Earth model.” Phys. Earth Plan. Int. 25:297-356.

*  http://ds.iris.edu/ds/products/emc-prem/

* https://www.cfa.harvard.edu/~lzeng/papers/PREM.pdf


To use the PREM model, you should add the following to your input file:

.. code:: yaml

  material:
    kind: prem


Note that the PREM model only makes sense when you simulate the Earth and therefore use the appropriate
axisymmetric domain for the Earth.