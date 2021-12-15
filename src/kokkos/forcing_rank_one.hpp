/*
//@HEADER
// ************************************************************************
//
// rank_one_forcing.hpp
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

#ifndef RANK_ONE_FORCING_HPP_
#define RANK_ONE_FORCING_HPP_

#include "KokkosBlas1_fill.hpp"
#include "KokkosBlas1_scal.hpp"

template <typename sc_t, typename state_d_t>
class RankOneForcing
{
  static_assert
  (is_kokkos_1dview< state_d_t >::value,
   "Rank-1 forcing needs a rank-1 kokkos view state");

  /*
    * on host, we precompute and store the full time signal
    since this is a just an array with as many elements as the number
    of time steps because the forcing source acts at a single grid point

    * on device we have an array with as many entries as number
    of velocity grid points and when we need to evaluate the forcing,
    we just copy from host a single value to the right location on device
   */

  using state_h_t = typename state_d_t::host_mirror_type;

  // f_h_ contains the full time series of the signal
  state_h_t f_h_;

  // f_d_ contains the forcing vector over the mesh
  state_d_t f_d_;

  // myVpGid_ identifies which velocity point the signal is located at
  // remember that forcing always acts on a velocity point, not a stress point.
  std::size_t myVpGid_ = 0;

  // max freq of all sources
  sc_t maxFreq_ = {};
  const sc_t dt_ = {};
  const std::size_t NSteps_ = {};

public:
  template <typename signal_t, typename parser_t, typename mesh_info_t, typename app_t>
  RankOneForcing(const signal_t & signalObj,
		 const parser_t & parser,
		 const mesh_info_t & meshInfo,
		 const app_t & appObj,
		 const sc_t depthKm,
		 const sc_t angleDeg)
    : f_h_("Fh", parser.getNumSteps()),
      f_d_("Fd", meshInfo.getNumVpPts()),
      dt_(parser.getTimeStepSize()),
      NSteps_(parser.getNumSteps()),
      maxFreq_(signalObj.getFrequency())
  {
    const auto gidsVp = appObj.viewGidListHost(dofId::vp);
    const auto coords = appObj.viewCoordsHost(dofId::vp);

    // find the vpGid identifying the grid point where the source is mapped to
    const auto domainSurfaceRadiusKm = meshInfo.getMaxRadiusKm();
    const sc_t myRadiusKm = domainSurfaceRadiusKm - depthKm;
    mapPointSourceToGridPoint(angleDeg, myRadiusKm, depthKm,
    			      meshInfo.viewDomainBounds(),
    			      meshInfo.getNumVpPts(), gidsVp, coords,
    			      meshInfo.getAngularSpacing(), myVpGid_);

    KokkosBlas::fill(f_h_, constants<sc_t>::zero());
    KokkosBlas::fill(f_d_, constants<sc_t>::zero());

    // compute and store the full time series of the signal on host
    storeSignalTimeSeries(signalObj);
  }

  template <typename parser_t, typename mesh_info_t, typename app_t>
  RankOneForcing(const parser_t    & parser,
		 const mesh_info_t & meshInfo,
		 const app_t	   & appObj)
    : RankOneForcing(Signal<sc_t>(parser.getSourceSignalKind(),
				  parser.viewDelays()[0],
				  parser.viewPeriods()[0]),
		     parser, meshInfo, appObj,
		     parser.viewDepths()[0],
		     parser.viewAngles()[0])
  {}

public:
  sc_t getMaxFreq() const{
    return maxFreq_;
  }

  std::size_t getVpGid() const{
    return myVpGid_;
  }

  sc_t getForcingValueAtStep(const std::size_t & step) const{
    return f_h_(step-1);
  }

  state_d_t viewForcingDevice() const{
    return f_d_;
  }

  void evaluate(const sc_t & time, const std::size_t & step)
  {
    const auto src = Kokkos::subview(f_h_, step-1);
    const auto des = Kokkos::subview(f_d_, myVpGid_);
    Kokkos::deep_copy(des, src);
  }

  // void complexityOfEvaluateMethod(double & memCostMB, double & flopsCost) const
  // {
  //   // // no operation is done during evaluate, just copying, see above
  //   const double memMBCostFill = 0.; //1.*( f_d_.extent(0)*sizeof(sc_t) )/1024./1024.;
  //   // for copy we have one read + one write
  //   const double memMBCostCopy = 1.*( 2*sizeof(sc_t) )/1024./1024.;
  //   memCostMB = memMBCostFill + memMBCostCopy;
  //   flopsCost = 0.;
  // }

private:
  template <typename signal_t>
  void storeSignalTimeSeries(const signal_t & signal)
  {
    // store the full time series of the signal into the host array
    sc_t time = constants<sc_t>::zero();
    for (std::size_t iStep = 1; iStep<=NSteps_; ++iStep)
    {
      signal(time, f_h_(iStep-1));
      time = iStep * dt_;
    }
  }
};

#endif
