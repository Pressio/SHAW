
# Creating and running a full case
Here we describe how to create and run a case. There are a few steps involved,
so we will describe each one.

Define an env var as follows:
```bash
export MYRUNDIR=${MYWORKDIR}/myFirstRun
mkdir ${MYRUNDIR}
```

## Generating the mesh
The code has been developed such that the mesh is generated with Python
and is used by the C++ code. There are two main reasons for this choice:
first, it allows us to decouple the mesh generation from the actual physics code;
second, we developed the code to support the concept of sample mesh,
which is a key feature for nonlinear ROMs. This is **not** needed right now
to solve the current elastic shear wave problem because this is a linear problem,
but it can be useful if, in the fugute, we extend the code to support nonlinear problems.

To specify the grid, one only needs to specify the grid for the velocity points because
the stress points are defined based on the staggered scheme (see paper).
Assume the test case you want to run uses a grid of 150 x 600 velocity points
along the radial and polar directions, respectively.
To generate the mesh files proceed as follows:
```python
cd ${ESWSRCDIR}/meshing
python create_single_mesh.py -nr 150 -nth 600 -working-dir ${MYRUNDIR}
```
This should generate inside `$MYRUNDIR}` a directory called `mesh150x600`
containing the following files;
```bash
-rw-r--r--  1 fnrizzi  staff   2.0M Aug 30 11:19 coeff_vp.dat
-rw-r--r--  1 fnrizzi  staff    12M Aug 30 11:19 graph_sp.dat
-rw-r--r--  1 fnrizzi  staff   7.2M Aug 30 11:19 graph_vp.dat
-rw-r--r--  1 fnrizzi  staff   229B Aug 30 11:19 mesh_info.dat
```
