/*
//@HEADER
// ************************************************************************
//
// material_setter.hpp
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

#ifndef MATERIAL_MAP_DISCONTINUITY_TO_MESH_HPP_
#define MATERIAL_MAP_DISCONTINUITY_TO_MESH_HPP_

template <typename mesh_info_t, typename gids_t, typename coords_t, typename labels_t>
void mapDiscontinuityToMesh(const mesh_info_t & meshInfo,
			    const gids_t  & gidsVp,
			    const gids_t  & gidsSp,
			    const coords_t & coordsVp,
			    const coords_t & coordsSp,
			    const labels_t & labelsSp,
			    const sc_t & requestedDiscontinuityDepthKm,
			    sc_t & mappedDiscontinuityDepthKm,
			    int_t & vpGidDiscont,
			    int_t & spGidDiscont)
{
  //     std::cout << "\nMapping discontinuity layer to grid" << std::endl;

  //     constexpr auto one	  = constants<sc_t>::one();
  //     constexpr auto thous  = constants<sc_t>::thousand();
  //     constexpr auto esrMeters = constants<sc_t>::earthSurfaceRadiusMeters();
  //     constexpr auto esrKm  = constants<sc_t>::earthSurfaceRadiusKm();

  //     // store the requested depth in Km of the discontinuity
  //     discontinuityRequestedDepthKm_ = parser_.depth2_;
  //     const auto discontinuityRequestedDepthM_ = parser_.depth2_ * thous;
  //     std::cout << "The requested depth = " << discontinuityRequestedDepthKm_ << " (km)" << std::endl;

  //     // find the first Vp GID of the grid point as close to the target depth as possible
  //     const auto numGptVp = meshInfo.getNumVpPts();
  //     auto trialDelta = std::numeric_limits<sc_t>::max() - discontinuityRequestedDepthM_;
  //     for (int_t iPt=0; iPt < numGptVp; ++iPt){
  //       const auto & ptGID      = gidsVp(iPt);
  //       const auto thisPtRadius = one/coordsVp(ptGID, 1); // meters
  //       const auto thisPtDepth  = esrMeters - thisPtRadius; //meters
  //       const auto delta	      = std::abs(thisPtDepth - discontinuityRequestedDepthM_);
  //       // here we need to check < only because we want to track the first Vp point
  //       // to meet the condition, otherwise with <= we would find the last point on
  //       // the curve along theta
  //       if ( delta < trialDelta ){
  // 	trialDelta = delta;
  // 	vpGidDiscont_ = ptGID;
  //       }
  //     }
  //     std::cout << " found vpPt " << vpGidDiscont_ << std::endl;

  //     // find the first srp GID of the grid point as close to the target depth as possible
  //     const auto numGptSp = meshInfo.getNumSpPts();
  //     trialDelta = std::numeric_limits<sc_t>::max() - discontinuityRequestedDepthM_;
  //     for (int_t iPt=0; iPt < numGptSp; ++iPt){
  //       const auto & ptGID      = gidsSp(iPt);
  //       const auto thisPtRadius = one/coordsSp(ptGID, 1); // meters
  //       const auto thisPtDepth  = esrMeters - thisPtRadius; //meters
  //       const auto delta      = std::abs(thisPtDepth - discontinuityRequestedDepthM_);
  //       if (labelsSp(ptGID)==2){
  // 	// here we need to check < only because we want to track the first point
  // 	// to meet the condition, otherwise with <= we would find the last point on
  // 	// the curve along theta
  // 	if ( delta < trialDelta ){
  // 	  trialDelta = delta;
  // 	  spGidDiscont_ = ptGID-1; // want the gid of the srp point, which is just one before
  // 	}
  //       }
  //     }
  //     std::cout << " found spPt " << spGidDiscont_ << std::endl;

  //     // the gid of the srp point can also be found using grid geometry as follows
  //     const int_t den = meshInfo.getNumPtsAlongTheta()-1;
  //     const int_t jumps = int_t(std::floor(vpGidDiscont_/den));
  //     const auto spGidDiscont2 = (jumps*den) + vpGidDiscont_;
  //     std::cout << " found spPt22 " << spGidDiscont2 << std::endl;
  //     if (spGidDiscont2 != spGidDiscont_){
  //       throw std::runtime_error("Something wrong when mapping discontinuity, \
  // the gids compute for srp do not match from the two methods used.");
  //     }

  //     // find the mapped discontinuity
  //     discontinuityMappedDepthKm_      = (esrMeters - (one/coordsVp(vpGidDiscont_,1)))/thous;
  //     discontinuityMappedRadiusKm_     = esrKm - discontinuityMappedDepthKm_;
  //     discontinuityMappedRadiusMeters_ = discontinuityMappedRadiusKm_ * thous;

  //     std::cout << "Discontinuity depth mapped to depth (Km) "
  // 	      << std::setprecision(dblFmt) << discontinuityMappedDepthKm_
  // 	      << " with radius (km) = " << discontinuityMappedRadiusKm_
  // 	      << " with VpGid = " << vpGidDiscont_
  // 	      << " with SpGid = " << spGidDiscont_
  // 	      << std::endl;
}

#endif
