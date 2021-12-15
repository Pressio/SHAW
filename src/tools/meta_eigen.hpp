/*
//@HEADER
// ************************************************************************
//
// meta_eigen.hpp
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

#ifndef METAF_EIGEN_HPP_
#define METAF_EIGEN_HPP_

#include <type_traits>
#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/SparseCore"

template <typename T, typename enable = void>
struct is_vector_eigen : std::false_type {};

template <typename T>
struct is_vector_eigen<
  T,
  typename std::enable_if<
    std::is_same<
      T, Eigen::Matrix<typename T::Scalar, Eigen::Dynamic, 1>
      >::value
    >::type
  > : std::true_type{};


template <typename T, typename enable = void>
struct is_col_major_dynamic_matrix_eigen : std::false_type {};

template <typename T>
struct is_col_major_dynamic_matrix_eigen<
  T,
  typename std::enable_if<
    std::is_same<
      T, Eigen::Matrix<typename T::Scalar, Eigen::Dynamic, Eigen::Dynamic>
      >::value
    >::type
  > : std::true_type{};


template <typename T, typename enable = void>
struct is_row_major_dynamic_matrix_eigen : std::false_type {};

template <typename T>
struct is_row_major_dynamic_matrix_eigen<
  T,
  typename std::enable_if<
    std::is_same<
      T,
      Eigen::Matrix<typename T::Scalar,
		    Eigen::Dynamic, Eigen::Dynamic,
		    Eigen::RowMajor>
      >::value
    >::type
  > : std::true_type{};


template <typename T, typename enable = void>
struct is_dynamic_matrix_eigen : std::false_type {};

template <typename T>
struct is_dynamic_matrix_eigen<
  T,
  typename std::enable_if<
    is_row_major_dynamic_matrix_eigen<T>::value or
    is_col_major_dynamic_matrix_eigen<T>::value
    >::type
  > : std::true_type{};

#endif
