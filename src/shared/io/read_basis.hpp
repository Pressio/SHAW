/*
//@HEADER
// ************************************************************************
//
// read_basis.hpp
//                     		Pressio/SHAW
//                         Copyright 2019
// National Technology & Engineering Solutions of Sandia, LLC (NTESS)
//
// Under the terms of Contract DE-NA0003525 with NTESS, the
// U.S. Government retains certain rights in this software.
//
// Pressio is licensed under BSD-3-Clause terms of use:
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
// contributors may be used to endorse or promote products derived
// from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Francesco Rizzi (fnrizzi@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef READ_BASIS_HPP_
#define READ_BASIS_HPP_

#include "matrix_read.hpp"

#ifdef SHW_HAVE_TPL_EIGEN
template <typename sc_t, typename int_t, typename dmat_t>
typename std::enable_if< is_dynamic_matrix_eigen<dmat_t>::value, dmat_t >::type
readBasis(const std::string fileName,
	  const int_t targetRomSize,
	  const int_t useBinary)
{
  /* here I might want only targetRomSize of basis vectors,
   * so to do this I have to:
   * 1. read the full basis vectors
   * 2. only extract the target columns I want
  */

  dmat_t M;
  if (useBinary == 1){
    readBinaryMatrixWithSize(fileName, M);
  }
  else{
    readAsciiMatrixWithSize(fileName, M);
  }

  // use the native functionalities to extract the target set of columns but
  // we cannot just use the native functionalities since Eigen uses expressions
  // to represnet things. So we need to construct a new object and return it.
  return dmat_t( M.leftCols(targetRomSize) );
}
#endif

template <typename sc_t, typename int_t, typename dmat_t>
typename std::enable_if< is_col_major_matrix_kokkos<dmat_t>::value, dmat_t >::type
readBasis(const std::string fileName,
	  const int_t targetRomSize,
	  const int_t useBinary)
{
  static_assert( has_host_space<dmat_t>::value,
		 "readBasis: you need to read on host memory first");

  // /* here I might want only targetRomSize of basis vectors,
  //  * so to do this I have to:
  //  * 1. read the full basis vectors
  //  * 2. only extract the target columns I want
  // */
  dmat_t M("M",1,1);
  if (useBinary == 1){
    readBinaryMatrixWithSize(fileName, M);
  }
  else{
    readAsciiMatrixWithSize(fileName, M);
  }

  // only return target subset of columns
  std::pair<std::size_t, std::size_t> indPair(0, targetRomSize);
  return Kokkos::subview(M, Kokkos::ALL, indPair);
}

#endif
