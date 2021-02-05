
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
