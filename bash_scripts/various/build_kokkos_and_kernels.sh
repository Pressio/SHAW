#!/usr/bin/env bash

set -e

function buildKokkos(){
    local kind=$1
    local arch=$2
    local ompBool=Off
    [[ ${kind} == serial ]] && ompBool=Off || ompBool=On

    [[ ! -d ${KOKKOS_BUILD_DIR} ]] && mkdir -p ${KOKKOS_BUILD_DIR}
    cd ${KOKKOS_BUILD_DIR} #&& rm -rf CMakeCache* core/*

    if [[ ${arch} == none ]]; then
	cmake -DCMAKE_CXX_COMPILER=${CXX} \
	      -DCMAKE_BUILD_TYPE="Release" \
	      -DCMAKE_INSTALL_PREFIX=${KOKKOSPFX} \
	      -DKokkos_ENABLE_TESTS=Off \
	      -DKokkos_ENABLE_OPENMP=${ompBool} \
	      -DKokkos_ENABLE_AGGRESSIVE_VECTORIZATION=Off \
	      ${KOKKOS_SRC}
    else
	cmake -DCMAKE_CXX_COMPILER=${CXX} \
	      -DCMAKE_BUILD_TYPE="Release" \
	      -DCMAKE_INSTALL_PREFIX=${KOKKOSPFX} \
	      -DKokkos_ARCH_${arch}=On \
	      -DKokkos_ENABLE_TESTS=Off \
	      -DKokkos_ENABLE_OPENMP=${ompBool} \
	      -DKokkos_ENABLE_AGGRESSIVE_VECTORIZATION=Off \
	      ${KOKKOS_SRC}
    fi
    make -j4
    make install
    cd ..
}

function buildKokkosKernels(){
    echo "Kernels using the KokkosPFX= ${PFX_DIR}"
    local kind=$1
    local arch=$2

    echo $BLAS_ROOT
    if [[ -z ${BLAS_ROOT} ]]; then
	echo "BLAS_ROOT env var must be defined"
	exit 1
    else
	echo "Kernels using BLAS_ROOT=${BLAS_ROOT}"
    fi

    echo $BLASLIBNAME
    if [[ -z ${BLASLIBNAME} ]]; then
	echo "BLASLIBNAME env var must be defined"
	exit 1
    else
	echo "Kernels using BLASLIBNAME=${BLASLIBNAME}"
    fi

    echo $LAPACK_ROOT
    if [[ -z ${LAPACK_ROOT} ]]; then
	echo "LAPACK_ROOT env var must be defined"
	exit 1
    else
	echo "Kernels using LAPACK_ROOT=${LAPACK_ROOT}"
    fi

    echo $LAPACKLIBNAME
    if [[ -z ${LAPACKLIBNAME} ]]; then
	echo "LAPACKLIBNAME env var must be defined"
	exit 1
    else
	echo "Kernels using LAPACKLIBNAME=${LAPACKLIBNAME}"
    fi

    [[ ! -d ${KOKKOSKER_BUILD_DIR} ]] && mkdir -p ${KOKKOSKER_BUILD_DIR}
    cd ${KOKKOSKER_BUILD_DIR} && rm -rf CMakeCache* src/*

    cmake \
	-DCMAKE_VERBOSE_MAKEFILE=On \
	-DCMAKE_CXX_COMPILER=${CXX} \
	-DCMAKE_BUILD_TYPE="Release" \
	-DCMAKE_INSTALL_PREFIX=${KOKKOSKERPFX} \
	\
	-DKokkosKernels_ENABLE_TPL_LAPACK=On \
	-DKokkosKernels_LAPACK_ROOT=${LAPACK_ROOT} \
	-DLAPACK_LIBRARIES="${LAPACKLIBNAME}" \
	\
	-DKokkosKernels_ENABLE_TPL_BLAS=On \
	-DKokkosKernels_BLAS_ROOT=${BLAS_ROOT} \
	-DBLAS_LIBRARIES="${BLASLIBNAME}" \
	\
	-DKokkosKernels_INST_DOUBLE=On \
	-DKokkosKernels_INST_LAYOUTRIGHT=On \
	-DKokkosKernels_INST_LAYOUTLEFT=On \
	-DKokkosKernels_INST_ORDINAL_INT=Off \
	-DKokkosKernels_INST_ORDINAL_INT64_T=On \
	-DKokkosKernels_INST_OFFSET_INT=Off \
	-DKokkosKernels_INST_OFFSET_SIZE_T=On \
	\
	-DKokkosKernels_ENABLE_TESTS=Off \
	-DKokkos_ROOT=${KOKKOSPFX} \
	${KOKKOS_KER_SRC}

    make -j16
    make install
    cd ..
}

MYPWD=`pwd`

# create dir
[[ ! -d kokkos ]] && mkdir kokkos
cd kokkos

if [ ! -f kokkos-3.1.01.tar.gz ]; then
    cp ${ESWSRCDIR}/tpls/kokkos-3.1.01.tar.gz .
fi
if [ ! -d kokkos-3.1.01 ]; then
    tar zxf kokkos-3.1.01.tar.gz
fi

if [ ! -f kokkos-kernels-3.1.01.tar.gz ]; then
    cp ${ESWSRCDIR}/tpls/kokkos-kernels-3.1.01.tar.gz .
fi
if [ ! -d kokkos-kernels-3.1.01 ]; then
    tar zxf kokkos-kernels-3.1.01.tar.gz
fi

KOKKOS_SRC=${MYPWD}/kokkos/kokkos-3.1.01
KOKKOS_KER_SRC=${MYPWD}/kokkos/kokkos-kernels-3.1.01
KOKKOS_BUILD_DIR=${MYPWD}/kokkos/kokkos_build
KOKKOSKER_BUILD_DIR=${MYPWD}/kokkos/kokkos_kernels_build

buildKokkos omp none
buildKokkosKernels omp none

cd ${MYPWD}
