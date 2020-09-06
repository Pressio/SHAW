
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
OMP_NUM_THREADS=4; OMP_PLACES=threads; OMP_PROC_BIND=spread; ./shwave_fom input.yaml
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
This will take about two minutes because currently the SVD
is not done very efficiently. The efficiency of this will be improved in a future version.
This SVD step should generate inside `${MYRUNDIR}` the following files:
```
lsv_vp    : left singular vectors for the velocity
sva_vp    : singular values for the velocity
lsv_sp    : left singular vectors for the stresses
sva_sp    : singular values for the stresses
```

## Run the ROM
Before running the ROM, we need to figure out the number of modes to use.
To this end, we use the cumulative energy of the singular values computed above
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

These values have already been set inside the `rom_input.yaml`.
Of course, you can play with different energies and set the modes in the `rom_input.yaml`.
So we can now run the ROM by doing:
```
cd ${MYRUNDIR}
cp ${ESWSRCDIR}/tutorialRunFiles/rank1_rom_run/rom_input.yaml ${MYRUNDIR}
ln -s ${MYWORKDIR}/build/shwave_rom .
OMP_NUM_THREADS=4; OMP_PLACES=cores; OMP_PROC_BIND=true; ./shwave_rom rom_input.yaml
mkdir ./rom; mv snaps_* rom
```
This should generate inside `${MYRUNDIR}/rom` the following files:
```bash
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
./reconstructFomState --podmodes=./lsv_vp_0 binary --romsize=191 \
	--romsnaps=./rom/snaps_vp binary --fsize=1 --outformat=ascii \
	--timesteps=7200 --samplingfreq=12 --outfileappend=vp
```
which should create a file `fomReconstructedState_timestep_7200_vp` with the reconstructed wavefield.
Make sure that above you use the correct arguments otherwise it will not work!
For example, if you change the number of modes you need to ensure you pass the correct arg.

Let's then plot the FOM, ROM and error fields:
```bash
cp ${ESWSRCDIR}/tutorialRunFiles/rank1_rom_rom/*.py ${MYRUNDIR}
python plotWavefield.py
```
which should generate three plot as follows:<br>
<img src="https://github.com/fnrizzi/ElasticShearWaves/blob/master/tutorialRunFiles/rank1_rom_run/wavefield_fom.png" width="33%">
<img src="https://github.com/fnrizzi/ElasticShearWaves/blob/master/tutorialRunFiles/rank1_rom_run/wavefield_rom.png" width="33%">
<img src="https://github.com/fnrizzi/ElasticShearWaves/blob/master/tutorialRunFiles/rank1_rom_run/wavefield_error.png" width="33%">
