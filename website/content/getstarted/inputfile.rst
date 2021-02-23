Input File
##########

:breadcrumb: {filename}/inputfile.rst
:summary: Input File
:date: 2021-02-22 11:00

The input file is organized into sections: *general*, *io*, *source*, *material*.
Below we describe each one separately.

For a full working example of input file, look for example at this `demo <{filename}/demos/rank1fom.rst>`_:

#####################
`1. General Section`_
#####################

The general section is *mandatory*.
Contains settings for, e.g., mesh location, time stepping, etc

.. code:: yaml

  general:
    meshDir: something 		     # full path to mesh directory
    dt: 1.			     # time step size in seconds
    finalTime: 150.		     # final simulation time in seconds
    checkNumericalDispersion: true   # enable/disable check for numerical dispersion
    checkCfl: true                   # enable/disable CFL check


|
|


################
`2. IO Section`_
################

Defines if and how to collect data. Specifically, the code supports data
collection in two forms, namely a snapshot matrix and/or seismogram.
The snaptshot matrix stores the state over the full mesh sampled with a specific
frequency during the simulation time. The seismogram, instead, stores the velocity
signal sampled with a target frequency only at target receivers on the domain surface.
These receiver locations are set in the input file by providing their angles in degrees.


.. code:: yaml

  io:
   snapshotMatrix:
     binary: true             # set false if you want to print ascii, true for binary
     velocity:
       freq: 1                # every how many time steps to sample velocity field
       fileName: snaps_vp     # filename to save snapshots to
     stress:
       freq: 1		      # every how many time steps to sample stresses
       fileName: snaps_sp     # filename to save snapshots to

   seismogram:
     binary: false            # set false to write ascii file, true for binary
     freq: 10                 # every how many time steps to sample stresses
     receivers: [5, 10, ...]  # comma-separated angles (degrees) of all receiver locations
			      # on surface where to collect seismograms


Note that the IO section is *optional*, following these rules:

* if the whole section is missing in the yaml input, data collection is disabled and no output files are generated

* if you only want the `snapshotMatrix` matrix, just enable that node and omit the seismogram

* if you only want the `seismogram`, enable that node and omit the one for the snapshot matrix

|
|

############################
`3. Source/forcing Section`_
############################

The forcing section is *mandatory* to define the source signal and parameters.
As anticipated in the highlights in `main page <{filename}/index.rst>`_, the code supports three scenarios:

* single forcing term

* multi-forcing simulation solved by running sequential rank-1 solves

* multi-forcing simulation solved by running sequential rank-2 solves

`3.1 Single Forcing Run`_
-------------------------

For a standard run with just a single forcing term, you can set up
the corresponding node in the yaml input as:

.. code:: yaml

  source:
    signal:
      kind: ricker   # choices: sinusoid, gaussDeriv, ricker
      depth: 1100.0  # depth of the source in Km
      period: 40.    # period of the signal in seconds
      delay: 10.0    # delay in seconds


`3.2 Multi-forcing simulation using rank-1`_
--------------------------------------------
This is the case where you are interested in simulating multiple forcing terms
within the same simulation, but want to use the rank-1 formulation, i.e. discrete
state and forcing term are stored using 1D arrays. Therefore, the code solves
all the realizations sequentially.

For example, suppose that you want to explore the wave dynamics for a source
with fixed kind, period and delay, but for multiple source depths.
To this end, you can just set the depth yaml field to be a comma-separated
list of target depths in Kilometers.


.. code:: yaml

  source:
    signal:
     kind: ricker
     depth: [1100., 550., 650., ...] # km
     period: 40.	             # seconds
     delay: 10.0		     # seconds


`3.3 Multi-forcing simulation using rank-2`_
--------------------------------------------
This is the case where you are interested in simulating multiple forcing terms
within the same simulation, and want to use the rank-2 formulation, i.e. the discrete
state and forcing term are stored using 2D arrays, allowing to solve sets
of relizations simultaneously.

For example, suppose that you want to explore the  wave dynamics for a source
with fixed kind, period and delay, but for multiple source depths.
To this end, you can just set the depth yaml field to be a comma-separated
list of target depths in Kilometers and specify a forcingSize.
The forcingSize defines how many realizations are solved at once using the rank-2 formulation.


.. code:: yaml

  source:
    signal:
     kind: ricker
     depth: [1100.0, 550, ...] # km
     period: 40.	       # seconds
     delay: 10.0               # seconds
     forcingSize: 4	       # forcingSize>=2 enables rank-2 solution


|
|

############################
`4. Material Model Section`_
############################

| The material section is *mandatory*, in the input file you have to choose one.
| We currently support a single, two-layer, the PREM, or custom material model.


`3.1 Single Layer Material Model`_
----------------------------------

A single medium with no discontinuities.
You can provide coefficients to define a quadratic parametrization of the density and shear velocity profile.
For more details, e.g. meaning and units, see `this <{filename}/getstarted/materialmodels.rst>`_.

.. code:: yaml

  material:
    kind: unilayer
    layer: {density: [2000., 0., 0.], velocity: [5000., 0., 0.]}


`3.2 Two-layer Material Model`_
-------------------------------

Represents a material model with two layers, separated by a
single discontinuity as shown in the figure below.
Both the density and shear velocity only have radial dependence.
You can provide coefficients to define a quadratic parametrization
of the density and shear velocity profile.
For more details, e.g. meaning and units, see `this <{filename}/getstarted/materialmodels.rst>`_.

.. code:: yaml

  material:
    kind: bilayer
    layer1: {density: [2000., 0., 0.004], velocity: [5000., 1., 0.05.]}
    layer2: {depth: 556, density: [100., 0.05, 0.01], velocity: [5000., 0., 0.]}


`3.3 The PREM Material Model`_
------------------------------

The PREM model is a radial model representing the average Earth properties, and one of the most
commonly adoptedo ones. For more details, check the following references:

.. code:: yaml

  material:
    kind: prem


`3.4 Custom Material Model`_
----------------------------

For this, fill the code at top of "main_fom.cc" to setup your custom material model.

.. code:: yaml

  material:
    kind: custom
