#!/bin/bash

set -e

# load global variables
source ${PWD}/global_vars.sh

# parse cline arguments
source ${PWD}/cmd_line_options.sh

# check that all basic variables are set, otherwise leave
check_minimum_vars_set

echo ""
echo "--------------------------------------------"
echo " current setting is: "
echo ""
print_global_vars
echo ""
echo "--------------------------------------------"

# set env if not already set
if [[ ! -z ${SETENVscript} ]]; then
    echo "loading environment from ${SETENVscript}"
    source ${SETENVscript}
    echo "PATH = $PATH"
else
    echo "--with-env-script NOT set, so we assume env is set already"
fi


# # check that blas and lapack are set
# if [[ -z ${BLAS_ROOT} ]]; then
#     echo "error: BLAS_ROOT must be found in the environment, exiting"
#     exit 2
# fi

# create working dir if not existing
[[ ! -d ${WORKINGDIR} ]] && mkdir -p ${WORKINGDIR}

# wipe everything if set to 1
[[ $WIPEEXISTING == yes ]] && wipe_existing_data_in_target_dir

#---------------------------
#---------------------------
# go to working dir
cd ${WORKINGDIR}

# get eigen
if [ ! -d ${WORKINGDIR}/tpls/eigen ]; then
    EIGENVERSION=3.3.7
    EIGENUNPACKEDDIRNAME=eigen-${EIGENVERSION}

    mkdir -p ${WORKINGDIR}/tpls/eigen
    cd ${WORKINGDIR}/tpls/eigen
    cp ${TOPDIR}/../tpls/eigen-${EIGENVERSION}.tar.gz .
    tar zxf eigen-${EIGENVERSION}.tar.gz
    mv ${EIGENUNPACKEDDIRNAME} eigen
    cd ${WORKINGDIR}
fi

# get yaml-cpp
if [ ! -d ${WORKINGDIR}/tpls/yamlcpp ]; then

    VNAME=yaml-cpp-0.6.3
    TARNAME=${VNAME}.tar.gz
    UNPACKNAME=yaml-cpp-yaml-cpp-0.6.3
    LINK=https://github.com/jbeder/yaml-cpp/archive/${TARNAME}
    YAMLCPPPFX=${WORKINGDIR}/tpls/yamlcpp/install
    SRC=${WORKINGDIR}/tpls/yamlcpp/${UNPACKNAME}

    mkdir -p ${WORKINGDIR}/tpls/yamlcpp
    cd ${WORKINGDIR}/tpls/yamlcpp
    wget ${LINK}
    tar xf ${TARNAME}
    mkdir build && cd build
    cmake -DCMAKE_INSTALL_PREFIX=${YAMLCPPPFX} \
	  -DYAML_CPP_BUILD_TESTS=OFF \
	  -DCMAKE_CXX_COMPILER=${CXX} \
	  ${SRC}
    make -j6 install
    cd ${WORKINGDIR}
else
    YAMLCPPPFX="${WORKINGDIR}/tpls/yamlcpp/install"
fi

#----------------------
# build wave code
#----------------------
EIGENPATH="${WORKINGDIR}/tpls/eigen/eigen"

#USEOMP=OFF
#[[ ${WITHOPENMP} == yes ]] && USEOMP=ON

KOKKOSKERDIR=
if [[ $ARCH == mac ]]; then
    KOKKOSKERDIR=${KOKKOSKERPFX}/lib/cmake/KokkosKernels
else
    KOKKOSKERDIR=${KOKKOSKERPFX}/lib64/cmake/KokkosKernels
fi

bdirname=build
if [[ ! -d ${bdirname} ]]; then
    mkdir ${bdirname}
fi
# enter
cd ${bdirname} && rm -rf *
cmake -DCMAKE_CXX_COMPILER=${CXX} \
      -DCMAKE_VERBOSE_MAKEFILE:BOOL=TRUE \
      -DCMAKE_BUILD_TYPE=Release \
      \
      -DEIGEN_INCLUDE_DIR=${EIGENPATH} \
      \
      -DYAMLCPP_INCLUDE_DIR=${YAMLCPPPFX}/include \
      -DYAMLCPP_LIB_DIR=${YAMLCPPPFX}/lib \
      \
      -DBLAS_LIB_DIR=${BLAS_ROOT}/lib \
      -DHAVE_KOKKOS:BOOL=ON \
      -DKokkosKernels_DIR=${KOKKOSKERDIR} \
      ${CPPSRC}
make -j4
cd ..

#-DHAVE_OMP:BOOL=${USEOMP}

#----------------------
# go back where we started
cd ${TOPDIR}
