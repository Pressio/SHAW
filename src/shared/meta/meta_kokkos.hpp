
#ifndef METAF_KOKKOS_HPP_
#define METAF_KOKKOS_HPP_

#include <type_traits>
#include "Kokkos_Core.hpp"
#include "KokkosSparse_CrsMatrix.hpp"

template <typename T, typename enable = void>
struct is_host_space : std::false_type {};

#ifdef KOKKOS_ENABLE_SERIAL
template <> struct is_host_space<Kokkos::Serial> : std::true_type{};
#endif

#ifdef KOKKOS_ENABLE_OPENMP
template <> struct is_host_space<Kokkos::OpenMP> : std::true_type{};
#endif
//--------------------------------------------


template <typename T, typename enable = void>
struct has_host_space : std::false_type {};

template <typename T>
struct has_host_space<
  T,
  typename std::enable_if<
    Kokkos::is_view<T>::value and
    std::is_same<
      typename T::traits::memory_space, Kokkos::HostSpace
      >::value
    >::type
  > : std::true_type{};


template <typename T, typename enable = void>
struct is_accessible_on_host : std::false_type {};

template <typename T>
struct is_accessible_on_host<
  T,
  typename std::enable_if<
    Kokkos::is_view<T>::value and
    Kokkos::SpaceAccessibility<
      typename T::traits::execution_space, Kokkos::HostSpace
      >::accessible
    >::type
  > : std::true_type{};
//--------------------------------------------

template <typename T, typename enable = void>
struct is_kokkos_view : std::false_type {};

template <typename T>
struct is_kokkos_view<
  T, // kokkos vector if (1)view and (2) has rank=1
  typename std::enable_if<
    Kokkos::is_view<T>::value
    >::type
  > : std::true_type{};

template <typename T, typename enable = void>
struct is_kokkos_1dview : std::false_type {};

template <typename T>
struct is_kokkos_1dview<
  T, // kokkos vector if (1)view and (2) has rank=1
  typename std::enable_if<
    Kokkos::is_view<T>::value && T::traits::rank==1
    >::type
  > : std::true_type{};


template <typename T, typename enable = void>
struct is_kokkos_2dview : std::false_type {};

template <typename T>
struct is_kokkos_2dview<
  T, // kokkos vector if (1)view and (2) has rank=2
  typename std::enable_if<
    Kokkos::is_view<T>::value && T::traits::rank==2
    >::type
  > : std::true_type{};


template <typename T, typename enable = void>
struct is_kokkos_3dview : std::false_type {};

template <typename T>
struct is_kokkos_3dview<
  T, // kokkos vector if (1)view and (2) has rank=3
  typename std::enable_if<
    Kokkos::is_view<T>::value && T::traits::rank==3
    >::type
  > : std::true_type{};
//--------------------------------------------

template <typename T, typename enable = void>
struct is_col_major_matrix_kokkos : std::false_type {};

template <typename T>
struct is_col_major_matrix_kokkos<
  T,
  typename std::enable_if<
    Kokkos::is_view<T>::value and
    T::traits::rank==2 and
    std::is_same<
      typename T::traits::array_layout, Kokkos::LayoutLeft
      >::value
    >::type
  > : std::true_type{};

#endif
