Input File
##########

:breadcrumb: {filename}/inputfile.rst
:summary: Input File: description, organization and things to pay attention to
:date: 2021-02-22 11:00

.. role:: yellow
    :class: m-text m-warning

.. role:: red
    :class: m-text m-danger

.. role:: blue
    :class: m-text m-info

.. role:: green
    :class: m-text m-success


The input file is organized into sections: *general*, *io*, *source*, *material*.
One input file must contain a single instance of each of these sections.
Below we describe each section separately.

For complete, valid examples of input files, look at `the full template file <https://github.com/fnrizzi/SHAW/blob/master/sampleInputFiles/template.yaml>`_, or the `demos' input files <https://github.com/fnrizzi/SHAW/blob/master/demos>`_.

#####################
`1. General Section`_
#####################

The section contains the following self-explanatory fields:

.. code:: yaml

  general:
    meshDir: something 		     # full path to mesh directory
    dt: 1.			     # time step size in seconds
    finalTime: 150.		     # final simulation time in seconds
    checkNumericalDispersion: true   # enable/disable check for numerical dispersion
    checkCfl: true                   # enable/disable CFL check

.. note-danger:: The general section is *mandatory*: do not forget it when you create the input file!

|

################
`2. IO Section`_
################

Defines if and how to collect data. Specifically, the code supports data
collection in two forms, namely a snapshot matrix and/or seismogram.
The snaptshot matrix stores the state over the full mesh sampled with a specific
frequency during the simulation. The seismogram, instead, stores the velocity
sampled with a target frequency only at specific receivers located on the domain surface.
The receivers' locations are set in the input file by providing their angles in degrees.
Receivers mimic the role of, e.g., real recording startions on the Earth surface.

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

.. note-success:: The IO section is *optional*:

  * if the full section is omitted in the input, data collection is disabled and no output files are generated

  * if you only want the `snapshotMatrix` matrix, just enable that node and omit the seismogram node

  * if you only want the `seismogram`, enable that node and omit the one for the snapshotMatrix node

|

############################
`3. Source/forcing Section`_
############################

Contains the parametrization of the source signal.

.. note-danger:: The forcing section is *mandatory*: do not forget it when you create the input file!

		 You need to choose one of these options:

		 * single forcing term

		 * multi-forcing simulation solved by running sequential rank-1 solves

		 * multi-forcing simulation solved by running sequential rank-2 solves

Below we discuss these options in detail.

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

For a full example of this scenario, see `the first demo <{filename}/rank1fom.rst>`_.


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

If, instead of the depth, you want to sample the period, you can fix the depth
and just provide a list of signal period samples to solve for.
If you provide a list of samples for both the period and depth, then the code
will use a tensor-product to define all cases.
For example, if you specify 20 depths and 10 periods, the code
will thus solve 200 trajectories.

For a full example of this rank-1 multi-forcing scenario, see `the second demo <{filename}/rank1fommulti.rst>`_.

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
:yellow:`Note that the forcingSize must be a divisor of number of target samples.`
For example, if you specify 20 depths, the forcingSize must be a divisor of 20.
Another example, if you specify 20 depths and 10 periods, the total number of trajectories
to compute is 200, so forcingSize must be a divisor of 200.

.. code:: yaml

  source:
    signal:
     # ...
     # same fields/options shown in 3.2 above
     # ...
     forcingSize: 4	       # forcingSize>=2 enables rank-2 solution


|

############################
`4. Material Model Section`_
############################

Last but not least, we have the material model parametrization.
You need to choose one of the options below.

.. note-danger:: The material model section is *mandatory*: do not forget it when you create the input file!


`3.1 Single Layer Material Model`_
----------------------------------
A single medium with no discontinuities.
You can provide coefficients to define a quadratic parametrization of the density and shear velocity profile.
For more details, e.g. meaning and units, see `this <{filename}/materialmodels.rst>`_.

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
For more details, e.g. meaning and units, see `this <{filename}/materialmodels.rst>`_.

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
