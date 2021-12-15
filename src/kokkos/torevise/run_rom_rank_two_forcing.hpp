/*
//@HEADER
// ************************************************************************
//
// run_rom_rank_two_forcing.hpp
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

#ifndef LEAP_FROG_RUN_ROM_RANK_TWO_FORCING_HPP_
#define LEAP_FROG_RUN_ROM_RANK_TWO_FORCING_HPP_

#include <numeric>
#include <KokkosBlas1_axpby.hpp>
#include "KokkosBlas3_gemm.hpp"

namespace kokkosapp{

template <typename sc_t, typename int_t>
void complexityRankTwoForcing(const int_t nVp,
			      const int_t nSp,
			      const int_t fSize,
			      double & memCostMB,
			      double & flopsCost)
{
  using comp_t = Complexity<sc_t>;

  std::array<double, 3> mem = {};
  std::array<double, 3> flops = {};

  // xRomVp = xRomVp + dt * ( romJvp * xRomSp )
  comp_t::gemm_beta_one(nVp, nSp, fSize, mem[0], flops[0]);

  // xRomVp = xRomVp + dt * phiVp^T * f
  comp_t::axpby(nVp, fSize, mem[1], flops[1]);

  // xRomSp = xRomSp + dt * ( romJsp * xRomVp )
  comp_t::gemm_beta_one(nSp, nVp, fSize, mem[2], flops[2]);

  // accumulate
  memCostMB   = std::accumulate(mem.begin(), mem.end(),   0.);
  flopsCost = std::accumulate(flops.begin(), flops.end(), 0.);
}


// template <typename sc_t, typename state_d_t>
// struct Foo
// {
//   sc_t T_;
//   sc_t dt_;
//   state_d_t xRomVp_d_;
//   state_d_t phiVpRhoInv_;

//   Foo(sc_t T, sc_t dt, state_d_t xRomVp_d, state_d_t phiVpRhoInv)
//     : T_(T), dt_(dt), xRomVp_d_(xRomVp_d), phiVpRhoInv_(phiVpRhoInv)
//   {}

//   KOKKOS_INLINE_FUNCTION
//   void operator() (const std::size_t & j) const
//   {
//     const sc_t freq = 1./50.;
//     const sc_t delayTime = 120.;

//     const auto tDiff   = (T_-delayTime);
//     const auto tDiffSq = tDiff*tDiff;
//     const auto expTerm = std::exp( -freq*tDiffSq );
//     const auto result = -2.*tDiff*freq*expTerm;

//     for (std::size_t i=0; i<xRomVp_d_.extent(0); ++i){
//       xRomVp_d_(i, j) += phiVpRhoInv_(i,j) * result;
//     }
//   }
// };


template <
  typename step_t, typename sc_t, typename forcing_t,
  typename observer_t, typename rom_jac_d_t, typename state_d_t
  >
void runRomRankTwoForcing(const step_t & numSteps,
			  const sc_t dt,
			  const state_d_t & phiVpRhoInv,
			  forcing_t & forcingObj,
			  observer_t & observerObj,
			  const rom_jac_d_t romJvp_d,
			  const rom_jac_d_t romJsp_d,
			  state_d_t xRomVp_d,
			  state_d_t xRomSp_d)
{
  assert(xRomVp_d.extent(1) == xRomSp_d.extent(1));
  KokkosBlas::fill(xRomVp_d, constants<sc_t>::zero());
  KokkosBlas::fill(xRomSp_d, constants<sc_t>::zero());

  const char ct_N    = 'N';
  constexpr auto one = constants<sc_t>::one();

  // I need this because KokkosBlas::axpy (see below in loop)
  // was giving me a compile error, so I need to use axpby
  Kokkos::View<sc_t*> kvOnes("kv1", xRomVp_d.extent(1));
  KokkosBlas::fill(kvOnes, one);

  // for timings
  Kokkos::Timer timer;
  double dataCollectionTime = {};
  std::array<double, 3> perfTimes = {1e32,0.,0.}; //min, max, total

  //****** LOOP ******//
  const auto startTime = std::chrono::high_resolution_clock::now();
  sc_t timeVp	  = constants<sc_t>::zero();
  for (auto iStep = 1; iStep<=numSteps; ++iStep)
  {
    if (iStep % 2000 == 0) std::cout << "ROM-KOKKOS-B step = " << iStep << std::endl;

    auto f = forcingObj.getForcingAtStep(iStep);

    // ----------------
    // 1. do velocity
    timer.reset();
    // xRomVp = xRomVp + dt * ( romJvp * xRomSp )
    KokkosBlas::gemm(&ct_N, &ct_N, dt, romJvp_d, xRomSp_d, one, xRomVp_d);
    // xRomVp = xRomVp + phiVpTRhoInv * dt * f // Note: here f contains dt already
    // Foo<sc_t, state_d_t> fnc(timeVp, dt, xRomVp_d, phiVpRhoInv);
    // Kokkos::parallel_for(xRomVp_d.extent(1), fnc);
    KokkosBlas::axpby(f, phiVpRhoInv, kvOnes, xRomVp_d);
    const double ct1 = timer.seconds();
    timer.reset();
    observerObj.observe(dofId::vp, iStep, xRomVp_d);
    dataCollectionTime += timer.seconds();

    // update time
    timeVp = iStep * dt;

    // ----------------
    // 2. do stress
    timer.reset();
    // xRomSp = xRomSp + dt * ( romJsp * xRomVp )
    KokkosBlas::gemm(&ct_N, &ct_N, dt, romJsp_d, xRomVp_d, one, xRomSp_d);
    const double ct2 = timer.seconds();
    timer.reset();
    observerObj.observe(dofId::sp, iStep, xRomSp_d);
    dataCollectionTime += timer.seconds();

    // ----------------
    // 3. timing vars
    const double time = ct1+ct2;
    perfTimes[0] = std::min(perfTimes[0], time);
    perfTimes[1] = std::max(perfTimes[1], time);
    perfTimes[2] += time;
  }

  const auto finishTime = std::chrono::high_resolution_clock::now();
  const std::chrono::duration<double> elapsed = finishTime - startTime;
  std::cout << "\nloopTime = " << std::fixed << std::setprecision(10) << elapsed.count();
  std::cout << "\ndataCollectionTime = " << std::fixed << std::setprecision(10)
	    << dataCollectionTime << std::endl;

  const auto fSize = xRomVp_d.extent(1);
  double memCostMB, flopsCost = 0.;
  complexityRankTwoForcing<sc_t>(xRomVp_d.extent(0), xRomSp_d.extent(0), fSize,
   				 memCostMB, flopsCost);
  printPerf(numSteps, perfTimes, memCostMB, flopsCost);
}

}//end namespace kokkosapp
#endif
