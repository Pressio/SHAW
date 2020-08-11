
#ifndef READ_BASIS_HPP_
#define READ_BASIS_HPP_

#include "matrix_read.hpp"

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
