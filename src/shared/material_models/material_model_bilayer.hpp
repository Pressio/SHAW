/*
//@HEADER
// ************************************************************************
//
// material_model_bilayer.hpp
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

#ifndef MATERIAL_MODEL_BILAYER_HPP_
#define MATERIAL_MODEL_BILAYER_HPP_

#include "material_model_base.hpp"

template<typename scalar_t, typename parser_t>
class BilayerMaterialModel final : public MaterialModelBase<scalar_t>
{
  using profile_params_t = typename parser_t::profile_params_t;
  using discont_depth_t  = typename parser_t::discont_depth_t;
  const discont_depth_t  & discontDepthsKm_;
  const profile_params_t & densityParams_;
  const profile_params_t & velocityParams_;
  const scalar_t domainSurfaceRadiusMeters_ = {};

public:
  template<typename mesh_info_t>
  BilayerMaterialModel(const parser_t & parser,
		       const mesh_info_t & meshInfo)
    : discontDepthsKm_(parser.viewDiscontinuityDepthsKm()),
      densityParams_(parser.viewDensityParametrization()),
      velocityParams_(parser.viewVelocityParametrization()),
      domainSurfaceRadiusMeters_(meshInfo.getMaxRadius())
  {}

  void computeAt(const scalar_t & radiusFromCenterMeters,
		 const scalar_t & angleRadians,
		 scalar_t & density,
		 scalar_t & vs) const final
  {
    constexpr auto thous     = constants<scalar_t>::thousand();
    const auto thisPtDepthM = domainSurfaceRadiusMeters_ - radiusFromCenterMeters;

    if ( thisPtDepthM < discontDepthsKm_[1]*thous){
      /* point belongs to layer1 */
      const auto & dM = thisPtDepthM;

      auto & rhoC = densityParams_[0];
      auto & vsC  = velocityParams_[0];
      density = rhoC[0] + rhoC[1]*dM + rhoC[2]*dM*dM;
      vs      = vsC[0]  + vsC[1]*dM  + vsC[2]*dM*dM;
    }
    else{
      /* point belongs to layer2 */

      //compute depth from discontDepth_
      const auto dM = thisPtDepthM - (discontDepthsKm_[1]*thous);
      auto & rhoC = densityParams_[1];
      auto & vsC  = velocityParams_[1];
      density = rhoC[0] + rhoC[1]*dM + rhoC[2]*dM*dM;
      vs      = vsC[0]  + vsC[1]*dM  + vsC[2]*dM*dM;
    }
  }
};

#endif
