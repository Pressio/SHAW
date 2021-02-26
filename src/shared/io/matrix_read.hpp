/*
//@HEADER
// ************************************************************************
//
// matrix_read.hpp
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

#ifndef MATRIX_READ_HPP_
#define MATRIX_READ_HPP_

#include <fstream>
#include <iomanip>
#include <iostream>
#include <ostream>

template<class dmat_t>
typename std::enable_if< is_col_major_matrix_kokkos<dmat_t>::value >::type
readBinaryMatrixWithSize(const std::string filename, dmat_t & M)
{
  using int_t = std::size_t;
  using sc_t  = typename dmat_t::value_type;
  std::ifstream fin(filename, std::ios::in | std::ios::binary);
  fin.exceptions(std::ifstream::failbit | std::ifstream::badbit);

  int_t rows={};
  int_t cols={};
  fin.read((char*) (&rows),sizeof(int_t));
  fin.read((char*) (&cols),sizeof(int_t));
  const auto nBytes = rows*cols*sizeof(sc_t);
  Kokkos::resize(M, rows, cols);
  fin.read( (char *) M.data(), nBytes );
  if (!fin){
    std::cout << std::strerror(errno) << std::endl;
    throw std::runtime_error("ERROR READING binary file");
  }
  else
    std::cout << fin.gcount() << " bytes read\n";

  fin.close();
}

template<class dmat_t>
typename std::enable_if<is_kokkos_3dview<dmat_t>::value>::type
readBinaryMatrixWithSize(const std::string filename, dmat_t & M)
{
  using int_t = std::size_t;
  using sc_t  = typename dmat_t::value_type;
  std::ifstream fin(filename, std::ios::in | std::ios::binary);
  fin.exceptions(std::ifstream::failbit | std::ifstream::badbit);

  int_t s1={};
  int_t s2={};
  int_t s3={};
  fin.read((char*) (&s1),sizeof(int_t));
  fin.read((char*) (&s2),sizeof(int_t));
  fin.read((char*) (&s3),sizeof(int_t));
  const auto nBytes = s1*s2*s3*sizeof(sc_t);
  Kokkos::resize(M, s1, s2, s3);
  fin.read( (char *) M.data(), nBytes );

  if (!fin){
    std::cout << std::strerror(errno) << std::endl;
    throw std::runtime_error("ERROR READING binary file");
  }
  else
    std::cout << fin.gcount() << " bytes read\n";

  fin.close();
}

template <typename dmat_t>
typename std::enable_if< is_kokkos_2dview<dmat_t>::value >::type
readAsciiMatrixWithSize(const std::string fileName, dmat_t M)
{
  static_assert( is_accessible_on_host<dmat_t>::value,
		 "readAsciiMatrixWithSize: the kokkos view must have HostSpace to read");

  std::ifstream source; source.open(fileName, std::ios_base::in);
  std::string line, colv;
  {
    std::getline(source, line);
    std::istringstream in(line);
    std::string col1, col2;
    in >> col1; in >> col2;
    Kokkos::resize(M, std::stoi(col1), std::stoi(col2) );
  }

  // then read the actual data
  std::size_t iRow = 0;
  while (std::getline(source, line) ){
    std::istringstream in(line);
    for (std::size_t j=0; j<M.extent(1); ++j){
      in >> colv;
      M(iRow, j) = atof(colv.c_str());
    }
    iRow++;
  }
  source.close();
}

#endif
