/*
//@HEADER
// ************************************************************************
//
// fom_complexities.hpp
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

#ifndef FOM_COMPLEXITIES_HPP_
#define FOM_COMPLEXITIES_HPP_

namespace kokkosapp{

/*
  rank-1 forcing
*/
template <typename sc_t, typename state_d_t, typename app_t, typename forcing_t>
typename std::enable_if< is_kokkos_1dview<state_d_t>::value >::type
complexityFom(const state_d_t xVp,
	      const state_d_t xSp,
	      const app_t & fomObj,
	      const forcing_t & forcingObj,
	      double & memCostMB,
	      double & flopsCost)
{
  using ord_t = typename app_t::jacobian_ord_type;
  using comp_t = Complexity<sc_t>;

  const auto nVp = xVp.extent(0);
  const auto nSp = xSp.extent(0);

  std::array<double, 3> memMB = {};
  std::array<double, 3> flops = {};

  // // for now leave out forcing complexity
  // forcingObj.complexityOfEvaluateMethod(memMB[0], flops[0]);

  // spmv: xVp = xVp + dt * Jvp * xSp
  const auto nnz_j_vp = fomObj.getJacobianNNZ(dofId::vp);
  comp_t::template spmv<ord_t>(nnz_j_vp, nVp, memMB[0], flops[0]);

  // mult: xVp = xVp + dt * rhoInv * f
  // (note that we specify beta=1 case
  comp_t::mult_beta_one(nVp, memMB[1], flops[1]);

  // spmv: xSp = xSp + dt * Jsp * xVp
  const auto nnz_j_sp = fomObj.getJacobianNNZ(dofId::sp);
  comp_t::template spmv<ord_t>(nnz_j_sp, nSp, memMB[2], flops[2]);

  memCostMB = std::accumulate(memMB.begin(), memMB.end(), 0.);
  flopsCost = std::accumulate(flops.begin(), flops.end(), 0.);
}


/*
   rank-2 forcing
*/
template <typename sc_t, typename state_d_t, typename app_t, typename forcing_t>
typename std::enable_if< is_kokkos_2dview<state_d_t>::value >::type
complexityFom(const state_d_t xVp,
	      const state_d_t xSp,
	      const app_t & fomObj,
	      const forcing_t & forcingObj,
	      double & memCostMB,
	      double & flopsCost)
{
  using ord_t = typename app_t::jacobian_ord_type;
  using comp_t = Complexity<sc_t>;

  assert( xVp.extent(1) == xSp.extent(1) );
  const auto nVp    = xVp.extent(0);
  const auto nSp    = xSp.extent(0);
  const auto fSize = xVp.extent(1);

  std::array<double, 3> memMB = {};
  std::array<double, 3> flops = {};

  // // for now leave out forcing complexity
  // forcingObj.complexityOfEvaluateMethod(memMB[0], flops[0]);

  // xVp = xVp + dt * Jvp * xSp
  const auto nnz_j_vp = fomObj.getJacobianNNZ(dofId::vp);
  comp_t::template spmm<ord_t>(nnz_j_vp, nVp, fSize, memMB[0], flops[0]);

  // for rank-2 we do a parallel for over number of forcing realizations
  // each thread read/writes about 4 words and performs 3 flops
  // this is a rough approximation
  memMB[1] = 4.*fSize*sizeof(sc_t)/1024./1024.;
  flops[1] = 3.*fSize;

  // spmm: xSp = xSp + dt * Jsp * xVp
  const auto nnz_j_sp = fomObj.getJacobianNNZ(dofId::sp);
  comp_t::template spmm<ord_t>(nnz_j_sp, nSp, fSize, memMB[2], flops[2]);

  memCostMB = std::accumulate(memMB.begin(), memMB.end(), 0.);
  flopsCost = std::accumulate(flops.begin(), flops.end(), 0.);
}

}//end namespace kokkosapp
#endif
