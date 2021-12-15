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

template <typename dmat_t>
void readBasis(const std::string fileName,
	       dmat_t & phi,
	       const int useBinary)
{
  static_assert(is_col_major_matrix_kokkos<dmat_t>::value,
		"readBasis: view must be col major");
  static_assert( has_host_space<dmat_t>::value,
		 "readBasis: you need to read on host memory first");

  const bool fileContainsExtents = true;

  /* here I might want only targetRomSize of basis vectors,
   * so to do this I have to:
   * 1. read the full basis vectors
   * 2. only extract the target columns I want
  */
  dmat_t M("M",1,1);
  if (useBinary == 1){
    fillMatrixFromBinary(fileName, M, fileContainsExtents);
  }
  else{
    fillMatrixFromAscii(fileName, M, fileContainsExtents);
  }
  std::cout << " TONI " << M.extent(0) << " " << M.extent(1) << std::endl;

  if (phi.extent(0) != M.extent(0)){
    std::cout << " PHI resize " << std::endl;
    Kokkos::resize(phi, M.extent(0), phi.extent(1));
  }
  std::cout << " PHI " << phi.extent(0) << " " << phi.extent(1) << std::endl;

  Kokkos::parallel_for(phi.extent(0),
		       KOKKOS_LAMBDA(const std::size_t i){
			 for (std::size_t j=0; j<phi.extent(1); ++j){
			   phi(i,j) = M(i,j);
			 }
		       });

  // // only return target subset of columns
  // std::pair<std::size_t, std::size_t> indPair(0, targetRomSize);
  // return Kokkos::subview(M, Kokkos::ALL, indPair);
}

template<class dmat_t>
typename std::enable_if< is_col_major_matrix_kokkos<dmat_t>::value >::type
readBasis(const std::string filename,
	  dmat_t & M,
	  const int useBinary,
	  bool fileContainsExtents)
{
  if (useBinary == 1){
    fillMatrixFromBinary(filename, M, fileContainsExtents);
  }
  else{
    fillMatrixFromAscii(filename, M, fileContainsExtents);
  }
}

#endif
