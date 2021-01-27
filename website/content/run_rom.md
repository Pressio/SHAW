
# Creating and running a rank-1 Galerkin ROM
Here we describe how to create and run a rank-1 Galerkin ROM simulation
in a reproductive scenario.
Specifically, this is a self-contained demo that shows the full steps:
first we run the FOM, compute the POD basis, then use the basis to execute the ROM.

## Preparing the env
Set the following env variables:
```bash
export ESWSRCDIR=<path-to-the-code-repository>
export MYWORKDIR=<the-same-work-directory-used-for-building-process>
export MYRUNDIR=${MYWORKDIR}/myRank1RomRun
mkdir ${MYRUNDIR}
```
Here we assume you already ran the [step-by-step guide](./docs/build.md)
to build the code and used `MYWORKDIR` as the working directory for that procedure,
so that `${MYWORKDIR}/build` contains all the executables.

## Generating the mesh
See the same section in [this tutorial](./run_fom.md) for the FOM for more details on the grid.
To generate the mesh files proceed as follows:
```python
cd ${ESWSRCDIR}/meshing
python create_single_mesh.py -nr 200 -nth 1000 -working-dir ${MYRUNDIR}
```
This should generate a directory `${MYRUNDIR}/mesh200x1000` containing:
```bash
-rw-r--r--  1 fnrizzi  staff   4.5M Aug 30 12:20 coeff_vp.dat
-rw-r--r--  1 fnrizzi  staff    28M Aug 30 12:20 graph_sp.dat
-rw-r--r--  1 fnrizzi  staff    16M Aug 30 12:20 graph_vp.dat
-rw-r--r--  1 fnrizzi  staff   231B Aug 30 12:20 mesh_info.dat
```

## Run the FOM
After creating the input file, we can now link the FOM executable and run:
```bash
cd ${MYRUNDIR}
cp ${ESWSRCDIR}/tutorialRunFiles/rank1_rom_run/fom_input.yaml ${MYRUNDIR}
ln -s ${MYWORKDIR}/build/shwave_fom .
OMP_NUM_THREADS=4 OMP_PLACES=threads OMP_PROC_BIND=spread ./shwave_fom fom_input.yaml
mkdir ./fom; mv coords_* seismogram_* snaps_* fom
```
This should generate inside `${MYRUNDIR}/fom` the following files:
```
coords_sp.txt : coordinates of the velocity grid points
coords_vp.txt : coordinates of the stresses grid points
seismogram_0  : seismogram at the receiver locations set in fom_input.yaml
snaps_vp_0    : snapshot matrix for the velocity
snaps_sp_0    : snapshot matrix for the stresses
```

## Compute POD modes
After running the FOM, use the snapshot matrices to compute POD modes:
```bash
cd ${MYRUNDIR}
ln -s ${MYWORKDIR}/build/computeThinSVD .
./computeThinSVD ./fom 1 0
```
This will take a couple minutes because currently the SVD
is not done very efficiently. The current version of the SVD focuses on accuracy.
The efficiency of this will be improved in a future version, but since this
step is orthogonal to the rest, one can easily swap it with something else.
This SVD step should generate inside `${MYRUNDIR}` the following files:
```
lsv_vp    : left singular vectors for the velocity
sva_vp    : singular values for the velocity
lsv_sp    : left singular vectors for the stresses
sva_sp    : singular values for the stresses
```

## Run the ROM
Before running the ROM, we need to set the number of modes to use.
To this end, we use the cumulative energy of the singular values computed via SVD above.
We proceed as follows:
```bash
cd ${MYRUNDIR}
cp $ESWSRCDIR/python_scripts/cumulative_energy.py .
python cumulative_energy.py -working-dir . -percent 99.9999999
```
this should print something like:
```bash
target 0.999999999
Nvp: 191 Nsp: 180
```
where it says that to cover 99.9999999 % of the energy for the velocity
we need 191 modes, while for the stresses we need 180 modes.

These values for the number of modes have already been set inside the `rom_input.yaml`.
Of course, you can play with different energies and change
the modes in the `rom_input.yaml` accordingly.
We can now run the ROM by doing:
```bash
cd ${MYRUNDIR}
cp ${ESWSRCDIR}/tutorialRunFiles/rank1_rom_run/rom_input.yaml ${MYRUNDIR}
ln -s ${MYWORKDIR}/build/shwave_rom .
OMP_NUM_THREADS=4 OMP_PLACES=cores OMP_PROC_BIND=true ./shwave_rom rom_input.yaml
mkdir ./rom; mv snaps_* rom
```
This should be a quite fast run, about a couple seconds total,
and should create inside `${MYRUNDIR}/rom` the following files:
```
snaps_vp_0    : snapshot matrix for the velocity generalized coordinates
snaps_sp_0    : snapshot matrix for the stresses generalized coordinates
```

## Post-process data
To post-process the data, we first extract the wavefield at the final time
from the FOM snapshot as:
```bash
cd ${MYRUNDIR}
ln -s ${MYWORKDIR}/build/extractStateFromSnaps .
./extractStateFromSnaps --snaps=./fom/snaps_vp_0 binary \
	--fsize=1 --outformat=ascii --timesteps=7200 --samplingfreq=12 --outfileappend=vp
```
which should create a file `state_timestep_7200_vp` with the FOM wavefield.

Then we use the ROM results to reconstruct the wavefield at the same time by doing:
```bash
cd ${MYRUNDIR}
ln -s ${MYWORKDIR}/build/reconstructFomState .
./reconstructFomState --podmodes=./lsv_vp binary --romsize=191 \
	--romsnaps=./rom/snaps_vp_0 binary --fsize=1 --outformat=ascii \
	--timesteps=7200 --samplingfreq=12 --outfileappend=vp
```
which should create `fomReconstructedState_timestep_7200_vp` with the reconstructed wavefield.
Make sure that above you use the correct arguments otherwise it will not work!
For example, if you change the number of modes you need to ensure you pass the correct arg.

Let's now compute the error between the FOM and the ROM-reconstructed wavefield as:
```bash
cp $ESWSRCDIR/python_scripts/compute_state_error.py .
python compute_state_error.py -dof-name vp \
	-fom-state ${PWD}/state_timestep_7200_vp -approx-state ${PWD}/fomReconstructedState_timestep_7200_vp
```
which should print:
```
 vp_fom_minmax    = -2.3482838537911068e-07 2.1150382886843842e-07
 vp_approx_minmax = -2.3482604603492044e-07 2.1151956606326872e-07
 vp_err_minmax    = -3.3787625883888276e-10 9.36860518631533e-10
 vp_err_abs_rel_ltwo_norms = 1.2584414985364181e-08 0.0029526114722209086
 vp_err_abs_rel_linf_norms = 9.36860518631533e-10 0.003989553976275272
```
where the second to the last line shows `vp_err_abs_rel_ltwo_norms` reporting
the absolute (1.2584414985364181e-08) and relative (0.0029526114722209086) L2 error.
The result indicates an excellent relative error of about 0.3%.

Let's then plot the FOM, ROM and error fields:
```bash
cp ${ESWSRCDIR}/tutorialRunFiles/rank1_rom_run/*.py ${MYRUNDIR}
python plotWavefield.py
```
which should generate three plots as follows:<br>
<img src="https://github.com/fnrizzi/ElasticShearWaves/blob/master/tutorialRunFiles/rank1_rom_run/wavefield_fom.png" width="33%">
<img src="https://github.com/fnrizzi/ElasticShearWaves/blob/master/tutorialRunFiles/rank1_rom_run/wavefield_rom.png" width="33%">
<img src="https://github.com/fnrizzi/ElasticShearWaves/blob/master/tutorialRunFiles/rank1_rom_run/wavefield_error.png" width="33%">
