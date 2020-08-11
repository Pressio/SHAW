
#ifndef METAF_EIGEN_HPP_
#define METAF_EIGEN_HPP_

#include <type_traits>
#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/SparseCore"

//--------------------------------------------
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
//--------------------------------------------


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
