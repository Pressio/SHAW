
Performance
===========

The following plot shows performance results obtained for the
:ref:`rank-2 formulation <rank2fom>` on a workstation
with two 18-core Intel(R) Xeon(R) Gold 6154 CPU @ 3.00 GHz,
each with a 24.75MB L3 cache and 125GB total memory.
We enable hyperthreading, thus supporting a maximum of 36 logical threads per CPU,
so a total of 72 threads. We use GCC-8.3.1 and rely on kokkos
and kokkos-kernels version 3.1.01.
We use Blis-0.7.0 as the kokkos-kernels’ backend for all dense operations.
We use the OpenMP backend for Kokkos.

|

.. figure:: ../img/fom_cpu_ave.png
   :align: center
   :width: 95%

   M represents how many trajectories we are computing simultaneously:
   when M=1, this what we refer to as :ref:`rank-1 formulation <rank1fom>`,
   while M>=2 corresponds to what we refer to as :ref:`rank-2 formulation <rank2fom>`;
   N is the *total* number of dofs (velocities plus stresses) for the problem.


todo: put link to script to run performance test
