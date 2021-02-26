/*
//@HEADER
// ************************************************************************
//
// map_nominal_location_to_velocity_grid_point.hpp
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

#ifndef MAP_NOMINAL_LOCATION_TO_VELOCITY_GRID_POINT_HPP_
#define MAP_NOMINAL_LOCATION_TO_VELOCITY_GRID_POINT_HPP_

#include <iostream>

template <
  typename scalar_t, typename bounds_t, typename int_t, typename gids_list_t, typename coords_t
  >
int_t mapNominalLocationToVelocityGridPoint(const scalar_t nominalAngleDeg,
					    const scalar_t nominalRadiusM,
					    const bounds_t & domainBounds,
					    const int_t numGptVp,
					    const gids_list_t & vpGidsList,
					    const coords_t & vpCoords)
{
  int_t result = 0;

  constexpr auto one	  = constants<scalar_t>::one();
  constexpr auto thousand = static_cast<scalar_t>(1000);

  const auto domainThLRad = degToRad(domainBounds[0]);
  const auto domainThRRad = degToRad(domainBounds[1]);

  // the nomial location
  auto nominalAngleRad = degToRad(nominalAngleDeg);
  std::cout << "Target coords (r [m], th [deg], th [rad]) = "
	    << nominalRadiusM << " "
	    << nominalAngleDeg    << " "
	    << nominalAngleRad    << std::endl;

  // check that source angle is within domain
  if ((nominalAngleDeg < domainBounds[0]) or (nominalAngleDeg > domainBounds[1])){
    throw std::runtime_error("Target angle is outside theta axis");
  }
  if ((nominalRadiusM < domainBounds[2]) or (nominalRadiusM > domainBounds[3])){
    throw std::runtime_error("Target angle is outside radius axis");
  }

  // find the velocity dof grid point that is closest to desired source
  auto trialDistance = std::numeric_limits<scalar_t>::max();
  for (int_t iPt=0; iPt < numGptVp; ++iPt){
    // make sure to loop over the graph and get the GID from it
    // because when we deal with sample mesh, iPt != ptGID
    // since one has more state points than actual residual points
    const auto & ptGID      = vpGidsList(iPt);
    const auto thisPtTh     = vpCoords(ptGID, 0);
    const auto thisPtRadius = one/vpCoords(ptGID, 1);

    // polar coords distance from test point to target
    auto thisD = computePolarDistance(thisPtRadius, thisPtTh, nominalRadiusM, nominalAngleRad);
    if ( thisD <= trialDistance ){
      trialDistance = thisD;
      result = ptGID;
    }
  }
  const auto srcTh = vpCoords(result,0);
  const auto srcRR = one/vpCoords(result,1);

  std::cout << "Point mapped to GID = "<< result
  	    << " with depth (km) = "
  	    << std::setprecision(dblFmt) << (domainBounds[3] - srcRR)/thousand
  	    << ", angle (rad) = "
  	    << std::setprecision(dblFmt) << (vpCoords(result,0))
  	    << ", angle (deg) = "
  	    << std::setprecision(dblFmt) << radToDeg(srcTh)
  	    << std::endl;

  return result;
}

#endif
