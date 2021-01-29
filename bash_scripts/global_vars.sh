#!/bin/bash

# top dir where this lives or is sourced
TOPDIR=${PWD}

# source
CPPSRC=${TOPDIR}/../src/cpp

# the working dir
WORKINGDIR=

# yes/no wipe existing data content of target directory
WIPEEXISTING=no

# env script
SETENVscript=

# location of kokkos if passed by user
KOKKOSPFX=
KOKKOSKERPFX=

# where yaml-cpp parser is
YAMLCPPPFX=

# var to detect the os type [linux or mac]
ARCH=
if [[ $OSTYPE == *"darwin"* ]]; then
    ARCH=mac
else
    ARCH=linux
fi

# if we want OpenMP enabled
WITHOPENMP=no

function wipe_existing_data_in_target_dir(){
    echo "Wiping ${CPPWORKINGDIR}/{data_*, build}"
    rm -rf ${WORKINGDIR}/build ${WORKINGDIR}/data_*
}

function print_global_vars(){
    echo "TOPDIR         = $TOPDIR"
    echo "CPPSRC         = $CPPSRC"
    echo "WORKINGDIR     = $WORKINGDIR"
    echo "WIPEEXISTING   = ${WIPEEXISTING}"
    echo "SETENVscript   = $SETENVscript"
    echo "KOKKOSPFX      = $KOKKOSPFX"
    echo "KOKKOSKERNPFX  = $KOKKOSKERPFX"
    echo "ARCH           = $ARCH"
    echo "OMP            = $WITHOPENMP"
}

function check_minimum_vars_set(){
    if [[ -z $WORKINGDIR ]]; then
	echo "--working-dir is empty, must be set: exiting"
	exit 1
    fi
}
