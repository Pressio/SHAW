#!/bin/bash

#---------------------------------
# trilinos generators
#---------------------------------
function tril_shaxipp(){
    trilinos_build_type
    trilinos_link_type
    trilinos_verbose_makefile_on
    trilinos_mpi_c_cxx_compilers
    trilinos_mpi_fortran_on
    trilinos_tests_off
    trilinos_examples_off
    trilinos_kokkos_serial
    trilinos_openblaslapack
    trilinos_packages_for_pressio
}

function tril_shaxipp_openmp(){
    trilinos_build_type
    trilinos_link_type
    trilinos_verbose_makefile_on
    trilinos_mpi_c_cxx_compilers
    trilinos_mpi_fortran_on
    trilinos_tests_off
    trilinos_examples_off
    trilinos_kokkos_omp
    trilinos_openblaslapack
    trilinos_packages_for_pressio
}
