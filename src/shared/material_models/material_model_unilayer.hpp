
#ifndef MATERIAL_MODEL_UNILAYER_HPP_
#define MATERIAL_MODEL_UNILAYER_HPP_

#include "material_model_base.hpp"

template<typename scalar_t, typename parser_t>
class UnilayerMaterialModel final : public MaterialModelBase<scalar_t>
{
  using profile_params_t = typename parser_t::profile_params_t;
  const profile_params_t & densityParams_;
  const profile_params_t & velocityParams_;
  const scalar_t domainSurfaceRadiusMeters_ = {};

public:
  template<typename mesh_info_t>
  UnilayerMaterialModel(const parser_t & parser,
			const mesh_info_t & meshInfo)
    : densityParams_(parser.viewDensityParametrization()),
      velocityParams_(parser.viewVelocityParametrization()),
      domainSurfaceRadiusMeters_(meshInfo.getMaxRadius())
  {}

  // evaluate density and shear velocity at target location
  void computeAt(const scalar_t & radiusFromCenterMeters,
		 const scalar_t & angleRadians,
		 scalar_t & density,
		 scalar_t & vs) const final
  {
    const auto dM = domainSurfaceRadiusMeters_ - radiusFromCenterMeters;

    auto & rhoC = densityParams_[0];
    auto & vsC  = velocityParams_[0];
    density = rhoC[0] + rhoC[1]*dM + rhoC[2]*dM*dM;
    vs	    = vsC[0] + vsC[1]*dM + vsC[2]*dM*dM;
  }
};

#endif
