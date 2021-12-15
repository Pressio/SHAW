/*
//@HEADER
// ************************************************************************
//
// vector_write.hpp
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

#ifndef VECTOR_WRITE_HPP_
#define VECTOR_WRITE_HPP_

#include <fstream>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <vector>

namespace impl
{

template<typename sc_t, typename size_t>
void write_vector_to_binary(const std::string filename,
			    const sc_t * A,
			    size_t n,
			    bool printSize = true)
{
  std::ofstream out(filename, std::ios::out | std::ios::binary | std::ios::trunc);
  if (printSize){
    out.write((char*) (&n), sizeof(size_t));
  }

  out.write((char*) A, n*sizeof(sc_t) );
  out.close();
}

template <typename mat_t, typename size_t>
void write_vector_to_ascii(const std::string fileName,
			   const mat_t & A,
			   size_t n,
			   bool printSize = true)
{
  std::ofstream file; file.open(fileName);
  if (printSize){
    file << n << std::endl;
  }

  for (size_t i=0; i<n; i++){
    file << std::setprecision(dblFmt) << A(i) << " \n";
  }
  file.close();
}

}// impl namespace

template <typename T>
typename std::enable_if< is_kokkos_1dview<T>::value >::type
writeToFile(const std::string fileName,
	    const T & v,
	    const bool useBinary,
	    bool printSize = true)
{
  static_assert(is_accessible_on_host<T>::value,
		"cannot call writeToFile for view not accessible on host");

  if (useBinary == 1){
    impl::write_vector_to_binary(fileName, v.data(), v.extent(0), printSize);
  }
  else{
    impl::write_vector_to_ascii(fileName, v, v.extent(0), printSize);
  }
}

#endif
