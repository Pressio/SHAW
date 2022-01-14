FROM ubuntu:focal
# Declaring build variables
# Timezone
ARG TZ=Europe/Warsaw
# CMake Version
ARG CMAKE_VERSION=3.18.6
# Compilers
ARG GNU_VER=9
ARG CC=gcc-$GNU_VER
ARG CXX=g++-$GNU_VER
ARG GFORTRAN=gfortran-$GNU_VER
ARG CC_PATH=/usr/bin/gcc
ARG CXX_PATH=/usr/bin/g++
# Options
ARG ENABLE_OPENMP=On
# TPLs
ARG KOKKOS_VER=3.5.00
ARG KOKKOSKERNELS_VER=3.5.00
ARG YAMLCPP_VER=0.7.0
ARG SRC=/src
ARG LIB=/tpl
ARG KOKKOS_SRC=$SRC/kokkos
ARG KOKKOSKERNELS_SRC=$SRC/kokkoskernels
ARG YAMLCPP_SRC=$SRC/yamlcpp
ARG KOKKOS_LIB=$LIB/kokkos
ARG KOKKOSKERNELS_LIB=$LIB/kokkoskernels
ARG YAMLCPP_LIB=$LIB/yamlcpp

# Setting timezone
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# Creating in and out directories
RUN mkdir /in
RUN mkdir /out
RUN mkdir -p $KOKKOS_SRC
RUN mkdir -p $KOKKOSKERNELS_SRC
RUN mkdir -p $YAMLCPP_SRC
RUN mkdir -p $KOKKOS_LIB
RUN mkdir -p $KOKKOSKERNELS_LIB
RUN mkdir -p $YAMLCPP_LIB

# System update and packages installation
RUN apt-get update && apt-get upgrade -y
# Installing Utilities
RUN apt-get install -y wget git make hwloc python3-pip python-is-python3
# Installing OpenMPI
RUN apt-get install -y openmpi-bin openmpi-doc
# Installing Libraries
RUN apt-get install -y libopenblas-dev liblapack-dev

# CMake installation
RUN wget -O cmake.sh https://github.com/Kitware/CMake/releases/download/v$CMAKE_VERSION/cmake-$CMAKE_VERSION-Linux-x86_64.sh
RUN sh cmake.sh --skip-license --exclude-subdir --prefix=/usr/local/
RUN rm cmake.sh

# Compilers installation
RUN apt-get install -y $CC $CXX $GFORTRAN
RUN update-alternatives --install $CC_PATH gcc /usr/bin/$CC 10
RUN update-alternatives --install $CXX_PATH g++ /usr/bin/$CXX 10
RUN update-alternatives --install /usr/bin/gfortran gfortran /usr/bin/$GFORTRAN 10
RUN update-alternatives --install /usr/bin/cc cc $CC_PATH 20
RUN update-alternatives --set cc $CC_PATH
RUN update-alternatives --install /usr/bin/c++ c++ $CXX_PATH 20
RUN update-alternatives --set c++ $CXX_PATH
RUN update-alternatives --install /usr/bin/fortrann fortrann /usr/bin/gfortran 20
RUN update-alternatives --set fortrann /usr/bin/gfortran

# Building TPLs
# Building Kokkos
WORKDIR $KOKKOS_SRC
RUN wget -O $KOKKOS_VER.tar.gz https://github.com/kokkos/kokkos/archive/refs/tags/$KOKKOS_VER.tar.gz
RUN tar -xzf $KOKKOS_VER.tar.gz
WORKDIR $KOKKOS_SRC/kokkos-$KOKKOS_VER
RUN mkdir build && cd build && cmake .. -DCMAKE_CXX_COMPILER=$CXX_PATH -DCMAKE_INSTALL_PREFIX=$KOKKOS_LIB -DKokkos_ENABLE_SERIAL=$ENABLE_OPENMP -DKokkos_ENABLE_OPENMP=On -DKokkos_ARCH_HSW=On -DKokkos_HWLOC_DIR=/usr/bin -DKokkos_ENABLE_AGGRESSIVE_VECTORIZATION=Off -DKokkos_ENABLE_TESTS=On && make install -j 8 && make test -j 8
# Building Kokkoskernels
WORKDIR $KOKKOSKERNELS_SRC
RUN wget -O $KOKKOSKERNELS_VER.tar.gz https://github.com/kokkos/kokkos-kernels/archive/refs/tags/$KOKKOSKERNELS_VER.tar.gz
RUN tar -xzf $KOKKOSKERNELS_VER.tar.gz
WORKDIR $KOKKOSKERNELS_SRC/kokkos-kernels-$KOKKOSKERNELS_VER
RUN mkdir build && cd build && cmake .. -DCMAKE_CXX_COMPILER=$CXX_PATH -DCMAKE_INSTALL_PREFIX=$KOKKOSKERNELS_LIB -DKokkos_ROOT=$KOKKOS_LIB/lib/cmake/Kokkos -DKokkosKernels_ENABLE_TPL_LAPACK=On -DKokkosKernels_ENABLE_TPL_BLAS=On -DKokkosKernels_INST_DOUBLE=On -DKokkosKernels_INST_LAYOUTRIGHT=On -DKokkosKernels_INST_LAYOUTLEFT=On -DKokkosKernels_INST_ORDINAL_INT=Off -DKokkosKernels_INST_ORDINAL_INT64_T=On -DKokkosKernels_INST_OFFSET_INT=Off -DKokkosKernels_INST_OFFSET_SIZE_T=On -DKokkosKernels_ENABLE_TESTS=On && make install -j 8 && make test -j 8
# Building YAMLCPP
WORKDIR $YAMLCPP_SRC
RUN wget -O $YAMLCPP_VER.tar.gz https://github.com/jbeder/yaml-cpp/archive/refs/tags/yaml-cpp-$YAMLCPP_VER.tar.gz
RUN tar -xzf $YAMLCPP_VER.tar.gz
WORKDIR $YAMLCPP_SRC/yaml-cpp-yaml-cpp-$YAMLCPP_VER
RUN mkdir build && cd build && cmake .. -DCMAKE_CXX_COMPILER=$CXX_PATH -DCMAKE_INSTALL_PREFIX=$YAMLCPP_LIB -DYAML_CPP_BUILD_TESTS=On && make install -j 8 && make test -j 8

# Clean up
WORKDIR /
RUN rm -rf $SRC

# Setting environment variables
ENV CC=$CC_PATH
ENV CXX=$CXX_PATH

# Setting workdir to /in
WORKDIR /in