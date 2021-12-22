#!/usr/bin/env bash

set -e

MYPWD=`pwd`

# where we need to work
TPLDIR=$1/tpls

kokkosversion=3.5.00
yamlcppversion=0.7.0

kokkosbackend=$2
arch=none


# --- nothing to addapt below --- #

function buildKokkos(){
    local arch=$1
    local backend=$2

    if [[ ${backend} == openmp ]]; then
	ompBool=On
    fi

    [[ ! -d ${KOKKOS_BUILD_DIR} ]] && mkdir -p ${KOKKOS_BUILD_DIR}
    cd ${KOKKOS_BUILD_DIR}

    if [[ ${arch} == none ]]; then
	cmake -DCMAKE_CXX_COMPILER=${CXX} \
	      -DCMAKE_BUILD_TYPE="Release" \
	      -DCMAKE_INSTALL_PREFIX=${KOKKOSPFX} \
	      -DKokkos_ENABLE_TESTS=Off \
	      -DKokkos_ENABLE_SERIAL=On \
	      -DKokkos_ENABLE_OPENMP=${ompBool} \
	      -DKokkos_ENABLE_AGGRESSIVE_VECTORIZATION=Off \
	      ${KOKKOS_SRC}
    else
	cmake -DCMAKE_CXX_COMPILER=${CXX} \
	      -DCMAKE_BUILD_TYPE="Release" \
	      -DCMAKE_INSTALL_PREFIX=${KOKKOSPFX} \
	      -DKokkos_ARCH_${arch}=On \
	      -DKokkos_ENABLE_TESTS=Off \
	      -DKokkos_ENABLE_SERIAL=On \
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
	-DKokkosKernels_ENABLE_TPL_BLAS=On \
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

function doKokkos
{
    cd ${TPLDIR}

    if [ ! -f kokkos-${kokkosversion}.tar.gz ]; then
	wget https://github.com/kokkos/kokkos/archive/refs/tags/${kokkosversion}.tar.gz
	mv ${kokkosversion}.tar.gz kokkos-${kokkosversion}.tar.gz
    fi
    if [ ! -d kokkos-${kokkosversion} ]; then
	tar zxf kokkos-${kokkosversion}.tar.gz
    fi

    if [ ! -f kokkos-kernels-${kokkosversion}.tar.gz ]; then
	wget https://github.com/kokkos/kokkos-kernels/archive/refs/tags/${kokkosversion}.tar.gz
	mv ${kokkosversion}.tar.gz kokkos-kernels-${kokkosversion}.tar.gz
    fi
    if [ ! -d kokkos-kernels-${kokkosversion} ]; then
	tar zxf kokkos-kernels-${kokkosversion}.tar.gz
    fi

    # sources for kokkos
    KOKKOS_SRC=${TPLDIR}/kokkos-${kokkosversion}
    KOKKOS_KER_SRC=${TPLDIR}/kokkos-kernels-${kokkosversion}
    # build dirs
    KOKKOS_BUILD_DIR=${TPLDIR}/kokkos-build
    KOKKOSKER_BUILD_DIR=${TPLDIR}/kokkos-kernels-build
    # prefixes
    KOKKOSPFX=${TPLDIR}/kokkos-install
    KOKKOSKERPFX=${TPLDIR}/kokkos-kernels-install

    buildKokkos ${arch} ${kokkosbackend}
    buildKokkosKernels ${arch} ${kokkosbackend}
}

function doYamlCpp()
{
    cd ${TPLDIR}

    VNAME=yaml-cpp-${yamlcppversion}
    TARNAME=${VNAME}.tar.gz
    UNPACKNAME=yaml-cpp-yaml-cpp-${yamlcppversion}
    LINK=https://github.com/jbeder/yaml-cpp/archive/${TARNAME}

    if [ ! -f ${TARNAME} ]; then
	wget ${LINK}
    fi
    if [ ! -d ${UNPACKNAME} ]; then
	tar zxf ${TARNAME}
    fi

    YAMLCPPPFX=${TPLDIR}/yamlcpp-install
    YAMLCPPBUILD=${TPLDIR}/yamlcpp-build

    if [ ! -d ${YAMLCPPBUILD} ]; then
	mkdir ${YAMLCPPBUILD} && cd ${YAMLCPPBUILD}
	cmake -DCMAKE_INSTALL_PREFIX=${YAMLCPPPFX} \
     	      -DYAML_CPP_BUILD_TESTS=OFF \
	      -DCMAKE_CXX_COMPILER=${CXX} \
	      -DYAML_BUILD_SHARED_LIBS=OFF \
	      ../${UNPACKNAME}
	make -j4 install
    fi
}


[[ ! -d ${TPLDIR} ]] && mkdir -p ${TPLDIR}
cd ${TPLDIR}

doKokkos
doYamlCpp

cd ${MYPWD}
