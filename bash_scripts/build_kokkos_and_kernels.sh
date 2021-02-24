#!/usr/bin/env bash

set -e

function buildKokkos(){
    local arch=$1
    local backend=$2

    if [[ ${backend} == openmp ]]; then
	ompBool=On
    fi
    if [[ ${backend} == serial ]]; then
	ompSerial=On
    fi

    [[ ! -d ${KOKKOS_BUILD_DIR} ]] && mkdir -p ${KOKKOS_BUILD_DIR}
    cd ${KOKKOS_BUILD_DIR} #&& rm -rf CMakeCache* core/*

    if [[ ${arch} == none ]]; then
	cmake -DCMAKE_CXX_COMPILER=${CXX} \
	      -DCMAKE_BUILD_TYPE="Release" \
	      -DCMAKE_INSTALL_PREFIX=${KOKKOSPFX} \
	      -DKokkos_ENABLE_TESTS=Off \
	      -DKokkos_ENABLE_SERIAL=${ompSerial} \
	      -DKokkos_ENABLE_OPENMP=${ompBool} \
	      -DKokkos_ENABLE_AGGRESSIVE_VECTORIZATION=Off \
	      ${KOKKOS_SRC}
    else
	cmake -DCMAKE_CXX_COMPILER=${CXX} \
	      -DCMAKE_BUILD_TYPE="Release" \
	      -DCMAKE_INSTALL_PREFIX=${KOKKOSPFX} \
	      -DKokkos_ARCH_${arch}=On \
	      -DKokkos_ENABLE_TESTS=Off \
	      -DKokkos_ENABLE_SERIAL=${ompSerial} \
	      -DKokkos_ENABLE_OPENMP=${ompBool} \
	      -DKokkos_ENABLE_AGGRESSIVE_VECTORIZATION=Off \
	      ${KOKKOS_SRC}
    fi
    make -j4
    make install
    cd ..
}

function buildKokkosKernels(){
    echo "Kernels using the KokkosPFX= ${KOKKOSPFX}"
    local arch=$1

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
	-DKokkosKernels_LAPACK_LIBRARIES="${LAPACKLIBNAME}" \
	\
	-DKokkosKernels_ENABLE_TPL_BLAS=On \
	-DKokkosKernels_BLAS_ROOT=${BLAS_ROOT} \
	-DKokkosKernels_BLAS_LIBRARIES="${BLASLIBNAME}" \
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
kokkosver=3.3.01

# create dir
[[ ! -d kokkos ]] && mkdir kokkos
cd kokkos

if [ ! -f kokkos-${kokkosver}.tar.gz ]; then
    cp ${ESWSRCDIR}/tpls/kokkos-${kokkosver}.tar.gz .
fi
if [ ! -d kokkos-${kokkosver} ]; then
    tar zxf kokkos-${kokkosver}.tar.gz
fi

if [ ! -f kokkos-kernels-${kokkosver}.tar.gz ]; then
    cp ${ESWSRCDIR}/tpls/kokkos-kernels-${kokkosver}.tar.gz .
fi
if [ ! -d kokkos-kernels-${kokkosver} ]; then
    tar zxf kokkos-kernels-${kokkosver}.tar.gz
fi

KOKKOS_SRC=${MYPWD}/kokkos/kokkos-${kokkosver}
KOKKOS_KER_SRC=${MYPWD}/kokkos/kokkos-kernels-${kokkosver}
KOKKOS_BUILD_DIR=${MYPWD}/kokkos/kokkos_build
KOKKOSKER_BUILD_DIR=${MYPWD}/kokkos/kokkos_kernels_build

buildKokkos none $1
buildKokkosKernels none $1

cd ${MYPWD}
