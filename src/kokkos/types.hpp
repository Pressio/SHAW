/*
//@HEADER
// ************************************************************************
//
// types.hpp
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

#ifndef SHAXIPP_KOKKOS_TYPES_HPP_
#define SHAXIPP_KOKKOS_TYPES_HPP_

#include "state_observer.hpp"
#include "seismogram.hpp"
#include "rank_one_forcing.hpp"
#include "rank_two_forcing.hpp"

namespace kokkosapp{

struct commonTypes
{
  using scalar_type = double;

  // compose parser with the various sections
  using p_gs_t  = ParserGeneralSection<scalar_type>;
  using p_io_t  = ParserIoSection<scalar_type>;
  using p_mm_t  = ParserMaterialModel<scalar_type>;
  using p_ss_t  = ParserForcingSection<scalar_type>;
  using parser_type  = InputParser<p_gs_t, p_io_t, p_mm_t, p_ss_t>;

  // mesh info class
  using mesh_info_type = MeshInfo<scalar_type>;

#ifdef KOKKOS_ENABLE_CUDA
  using device_mem_space = Kokkos::Cuda;
#else
  using device_mem_space = Kokkos::HostSpace;
#endif

  // jacobian is a sparse matrix
  using jacobian_ord_type = typename mesh_info_type::ordinal_type;
  using jacobian_d_type = KokkosSparse::CrsMatrix<scalar_type, jacobian_ord_type, device_mem_space>;

  // state observer: to collect state snapshots
  using observer_type = StateObserver<scalar_type>;

  // seismogram: to collect data at specific receivers on surface
  using seismogram_type  = Seismogram<scalar_type>;
};

struct rank1Types : commonTypes
{
  using typename commonTypes::scalar_type;
  using typename commonTypes::parser_type;
  using typename commonTypes::mesh_info_type;
  using typename commonTypes::jacobian_ord_type;
  using typename commonTypes::jacobian_d_type;
  using typename commonTypes::observer_type;
  using typename commonTypes::seismogram_type;
  using typename commonTypes::device_mem_space;

  // state is a rank-1 view
  using state_d_type = Kokkos::View<scalar_type*, device_mem_space>;
  using state_h_type = typename state_d_type::host_mirror_type;

  // forcing
  using forcing_type = RankOneForcing<scalar_type, state_d_type>;
};

struct rank2Types : commonTypes
{
  using typename commonTypes::scalar_type;
  using typename commonTypes::parser_type;
  using typename commonTypes::mesh_info_type;
  using typename commonTypes::jacobian_ord_type;
  using typename commonTypes::jacobian_d_type;
  using typename commonTypes::observer_type;
  using typename commonTypes::seismogram_type;
  using typename commonTypes::device_mem_space;

  // state is a rank-2 view
  using state_d_type = Kokkos::View<scalar_type**, device_mem_space>;
  using state_h_type = typename state_d_type::host_mirror_type;

  // forcing
  using signal_type = Signal<scalar_type>;
  using signal_instances_h_type = Kokkos::View<signal_type*, Kokkos::HostSpace>;
  using forcing_type = RankTwoForcing<scalar_type, signal_instances_h_type>;
};

}// end namespace
#endif
