/*
//@HEADER
// ************************************************************************
//
// vector_read.hpp
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

#ifndef VECTOR_READ_HPP_
#define VECTOR_READ_HPP_

#include <fstream>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <vector>

template <typename T>
typename std::enable_if< is_kokkos_1dview<T>::value >::type
readAsciiVectorWithSize(const std::string fileName, T v)
{
  static_assert( is_accessible_on_host<T>::value,
		 "readAsciiVectorWithSize: the kokkos view must have HostSpace to read");

  std::ifstream source; source.open(fileName, std::ios_base::in);
  std::string line, colv;
  {
    std::getline(source, line);
    std::istringstream in(line);
    std::string col1;
    in >> col1;
    Kokkos::resize(v, std::stoi(col1));
  }

  // then read the actual data
  std::size_t iRow = 0;
  while (std::getline(source, line) ){
    std::istringstream in(line);
    in >> colv;
    v(iRow) = atof(colv.c_str());
    iRow++;
  }
  source.close();
}

#endif
