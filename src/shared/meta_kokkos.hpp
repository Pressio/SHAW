/*
//@HEADER
// ************************************************************************
//
// meta_kokkos.hpp
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
