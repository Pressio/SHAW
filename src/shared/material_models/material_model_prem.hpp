/*
//@HEADER
// ************************************************************************
//
// material_model_prem.hpp
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

#ifndef MATERIAL_MODEL_PREM_HPP_
#define MATERIAL_MODEL_PREM_HPP_

#include "material_model_base.hpp"

template<typename scalar_t, typename parser_t>
class PremMaterialModel final : public MaterialModelBase<scalar_t>
{
  using profile_params_t = typename parser_t::profile_params_t;

  const scalar_t premEarthCmbKm_     = static_cast<scalar_t>(3480);
  const scalar_t premEarthCmbMeters_ = static_cast<scalar_t>(3480000);
  const scalar_t premEarthSurfaceKm_     = static_cast<scalar_t>(6371);
  const scalar_t premEarthSurfaceMeters_ = static_cast<scalar_t>(6371000);

public:
  template<typename mesh_info_t>
  PremMaterialModel(const parser_t & parser,
		    const mesh_info_t & meshInfo)
  {
    // for PREM, we need to make sure the domain loaded from mesh
    // has specific boundd since PREM is valid for Earth
    // with r_CMB = 3480 and r_earthSurface = 6371
    const auto loadedDomainSurfaceRadiusKm = meshInfo.getMaxRadiusKm();
    const auto loadedDomainMinimumRadiusKm = meshInfo.getMinRadiusKm();

    const auto diff1 = std::abs(loadedDomainSurfaceRadiusKm - premEarthSurfaceKm_);
    const auto diff2 = std::abs(loadedDomainMinimumRadiusKm - premEarthCmbKm_);
    if (diff1 > 1e-13 or diff2 > 1e-13)
    {
      throw std::runtime_error
	("To use the PREM model, your Earth domain must have r_cmb=3480 km \
and r_surface=6371 km. It seems you are using a domain/mesh that does NOT match this.");
    }
  }

  void computeAt(const scalar_t & radiusFromCenterMeters,
		 const scalar_t & angleRadians,
		 scalar_t & rho,
		 scalar_t & vs) const final
  {
    // If you use the Preliminary reference Earth model (PREM)
    // for your own research, please refer to
    // Dziewonski, A.M., and D.L. Anderson. 1981. “Preliminary reference Earth model.” Phys. Earth Plan. Int. 25:297-356.
    // IRIS DMC (2011), Data Services Products: EMC, A repository of Earth models, https://doi.org/10.17611/DP/EMC.1.
    // in your publication.

    // https://www.cfa.harvard.edu/~lzeng/papers/PREM.pdf

    constexpr auto thousand  = constants<scalar_t>::thousand();
    const auto rKm = radiusFromCenterMeters/thousand;
    const auto x   = rKm/premEarthSurfaceKm_;
    const auto xSq = x*x;
    const auto xCu = x*x*x;

    if(rKm >= 6356.0){
      // crust
      rho = 2.6;
      vs = 3.2;
    }
    else if(rKm >= 6346.6 and rKm < 6356.0){
      // crust
      rho = 2.9;
      vs = 3.9;
    }

    else if(rKm >= 6291.0 and rKm < 6346.6){
      rho = 2.691  + 0.6924*x;
      vs  = 2.1519 + 2.3481*x;
    }

    else if(rKm >= 6151.0 and rKm < 6291.0){
      rho = 2.691  + 0.6924*x;
      vs  = 2.1519 + 2.3481*x;
    }

    else if(rKm >= 5971.0 and rKm < 6151.0){
      // transition zone
      rho = 7.1089 - 3.8045*x;
      vs  = 8.9496 - 4.4597*x;
    }

    else if(rKm >= 5771.0 and rKm < 5971.0){
      // transition zone
      rho = 11.2494 - 8.0298*x;
      vs  = 22.3512 - 18.5856*x;
    }

    else if(rKm >= 5701.0 and rKm < 5771.0){
      // transition zone
      rho = 5.3197 - 1.4836*x;
      vs  = 9.9839 - 4.9324*x;
    }

    else if(rKm >= 5600.0 and rKm < 5701.0){
      // lower mantle part 3
      rho = 7.9565 - 6.4761*x  + 5.5283*xSq - 3.0807*xCu;
      vs = 22.3459 - 17.2473*x - 2.0834*xSq + 0.9783*xCu;
    }

    else if(rKm >= 3630.0 and rKm < 5600.0){
      // lower mantle part 2
      rho = 7.9565 - 6.4761*x  + 5.5283*xSq  - 3.0807*xCu;
      vs = 11.1671 - 13.7818*x + 17.4575*xSq - 9.2777*xCu;
    }

    else if(rKm >= 3480.0 and rKm < 3630.0){
      // lower mantle part 1
      rho = 7.9565 - 6.4761*x + 5.5283*xSq - 3.0807*xCu;
      vs  = 6.9254 + 1.4672*x - 2.0834*xSq + 0.9783*xCu;
    }

    else if(rKm >= 1221.5 and rKm < 3480.0){
      // outer core
      rho = 12.5815 - 1.2638*x - 3.6426*xSq - 5.5281*xCu;
      vs = 0.0;
    }
    else if(rKm < 1221.5){
      // inner core
      rho = 13.0885 - 8.8381*xSq;
      vs  = 3.6678  - 4.4475*xSq;
    }

    // convert vs from km/s to m/s
    vs *= thousand;
    // convert rho from g/cm^3 to kg/m^3
    rho *= thousand;
  }
};

#endif
