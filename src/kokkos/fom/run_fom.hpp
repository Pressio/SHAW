/*
//@HEADER
// ************************************************************************
//
// run_fom.hpp
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

#ifndef LEAP_FROG_RUN_FOM_HPP_
#define LEAP_FROG_RUN_FOM_HPP_

#include "fom_update_kernels.hpp"
#include "fom_complexities.hpp"

namespace kokkosapp{

template <
  typename step_t,
  typename sc_t,
  typename app_t,
  typename forcing_t,
  typename observer_t,
  typename seismo_t,
  typename state_d_t
  >
void runFom(const step_t & numSteps,
	    const sc_t dt,
	    const app_t & fomObj,
	    forcing_t & forcingObj,
	    observer_t & observerObj,
	    seismo_t & seismoObj,
	    state_d_t xVp_d,
	    state_d_t xSp_d)
{
  // zero states
  KokkosBlas::fill(xVp_d, constants<sc_t>::zero());
  KokkosBlas::fill(xSp_d, constants<sc_t>::zero());

  const auto jacVp_d     = fomObj.viewJacobianDevice(dofId::vp);
  const auto jacSp_d     = fomObj.viewJacobianDevice(dofId::sp);
  const auto rhoInvVp_d  = fomObj.viewInvDensityDevice(dofId::vp);

  // create mirrors of states if we need to collect data
  auto xVp_h = Kokkos::create_mirror_view(xVp_d);
  auto xSp_h = Kokkos::create_mirror_view(xSp_d);
  const auto snapshotsCollectionEnabled = observerObj.enabled();
  const auto seismogramEnabled = seismoObj.enabled();

  // to collec timings
  Kokkos::Timer timer;
  double dataCollectionTime = {};
  std::array<double, 3> perfTimes = {1e32,0.,0.}; //min, max, total

  //****** LOOP ******//
  const auto startTime  = std::chrono::high_resolution_clock::now();
  sc_t timeVp = constants<sc_t>::zero();
  for (auto iStep = 1; iStep<=numSteps; ++iStep)
  {
    if (iStep % 2000 == 0) std::cout << "Doing step = " << iStep << std::endl;

    // compute forcing for current time
    timer.reset();
    forcingObj.evaluate(timeVp, iStep);
    const double ct1 = timer.seconds();

    // ----------------
    // 1. do velocity
    timer.reset();
    updateVelocity(dt, xVp_d, xSp_d, jacVp_d, rhoInvVp_d, forcingObj);
    const double ct2 = timer.seconds();
    timer.reset();
    if (snapshotsCollectionEnabled or seismogramEnabled){
      // deep copy already fences so no need to explicitly fence
      Kokkos::deep_copy(xVp_h, xVp_d);
    }
    else{
      Kokkos::fence();
    }
    observerObj.observe(dofId::vp, iStep, xVp_h);
    seismoObj.storeVelocitySignalAtReceivers(iStep, xVp_h);
    dataCollectionTime += timer.seconds();

    // update time
    timeVp = iStep*dt;

    // ----------------
    // 2. do stress
    timer.reset();
    updateStress(dt, xSp_d, xVp_d, jacSp_d);
    const double ct3 = timer.seconds();
    timer.reset();
    if (snapshotsCollectionEnabled){
      Kokkos::deep_copy(xSp_h, xSp_d);
    }
    else{
      Kokkos::fence();
    }
    observerObj.observe(dofId::sp, iStep, xSp_h);
    dataCollectionTime += timer.seconds();

    // ----------------
    // 3. timing vars
    const double time = ct1+ct2+ct3;
    perfTimes[0] = std::min(perfTimes[0], time);
    perfTimes[1] = std::max(perfTimes[1], time);
    perfTimes[2] += time;
  }

  const auto finishTime = std::chrono::high_resolution_clock::now();
  const std::chrono::duration<double> elapsed = finishTime - startTime;
  std::cout << "\nloopTime = " << std::fixed << std::setprecision(10) << elapsed.count();
  std::cout << "\ndataCollectionTime = " << std::fixed << std::setprecision(10)
	    << dataCollectionTime << std::endl;

  // compute complexity and print
  double memCostMB, flopsCost = 0.;
  complexityFom<sc_t>(xVp_d, xSp_d, fomObj, forcingObj, memCostMB, flopsCost);
  printPerf(numSteps, perfTimes, memCostMB, flopsCost);
}

}//end namespace kokkosapp
#endif
