#!/usr/bin/env bash

function build_openblas(){
    local PWD=`pwd`
    local PARENTDIR=$PWD

    #-----------------------------------
    # create dir
    [[ ! -d openblas ]] && mkdir openblas
    cd openblas

    if [ ! -f v0.3.10.tar.gz ]; then
	#wget https://github.com/xianyi/OpenBLAS/archive/v0.3.10.tar.gz
	cp ${ESWSRCDIR}/tpls/v0.3.10.tar.gz .
    fi

    if [ ! -d OpenBLAS-0.3.10 ]; then
	tar zxf v0.3.10.tar.gz
    fi

    echo ${PWD}
    if [ ! -d install ]; then
	cd OpenBLAS-0.3.10

	BFName="${PWD}/../build.txt"
	IFName="${PWD}/../install.txt"

	echo "Building OpenBLAS"
	make BINARY=64 HOSTCC=$CC > ${BFName} 2>&1
	echo "Build output written to ${BFName}"
	if grep -q "OpenBLAS build complete" "${BFName}"; then
	    echo "OpenBLAS built successfull"
	else
	    echo "OpenBLAS built unsuccessfull"
	    exit 44
	fi

	echo "Installing OpenBLAS"
	make PREFIX=${PWD}/../install install > ${IFName} 2>&1
	echo "Install output written to ${IFName}"
    fi
}

build_openblas