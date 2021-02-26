/*
//@HEADER
// ************************************************************************
//
// map_point_source_to_velocity_grid_point.hpp
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

#ifndef MAP_POINT_SOURCE_TO_VELOCITY_GRID_POINT_HPP_
#define MAP_POINT_SOURCE_TO_VELOCITY_GRID_POINT_HPP_

#include <iostream>

template <
  typename scalar_t,
  typename bounds_t,
  typename int_t,
  typename gids_list_t,
  typename coords_t
  >
void mapPointSourceToGridPoint(const scalar_t srcAngleDeg,
			       const scalar_t srcRadiusKm,
			       const scalar_t srcDepthKm,
			       const bounds_t & domainBounds,
			       const int_t numGptVp,
			       const gids_list_t & gidsList,
			       const coords_t & coords,
			       const scalar_t dth,
			       int_t & pointGid)
{
  constexpr auto one	  = constants<scalar_t>::one();
  constexpr auto thousand = static_cast<scalar_t>(1000);

  const auto domainThLRad = degToRad(domainBounds[0]);
  const auto domainThRRad = degToRad(domainBounds[1]);

  // the location of the source desired by user
  auto srcAngleRad = degToRad(srcAngleDeg);
  std::cout << "Target Source coords (depth [km], r [km], th [deg], th [rad]) = "
	    << srcDepthKm  << " "
	    << srcRadiusKm << " "
	    << srcAngleDeg    << " "
	    << srcAngleRad    << std::endl;

  // check that source angle is within domain
  if ((srcAngleDeg < domainBounds[0]) or (srcAngleDeg > domainBounds[1]))
    throw std::runtime_error("Source angle is outside theta axis");
  if ((srcRadiusKm*thousand < domainBounds[2]) or (srcRadiusKm*thousand > domainBounds[3]))
    throw std::runtime_error("Source angle is outside radius axis");

  // check if we are on one of the edges, and if so, remap to
  // leftEdge + dtheta or rightEdge -dtheta
  // because on the left/right edges we cannot compute cotangent
  if ( essentiallyEqual(srcAngleDeg, domainBounds[0]) ){
    std::cout << "Source angle cannot be on left edge, mapping to: "
	      << std::setprecision(dblFmt) << domainThLRad+dth << std::endl;
    srcAngleRad = domainThLRad+dth;
  }
  if ( essentiallyEqual(srcAngleDeg, domainBounds[1]) ){
    std::cout << "Source angle cannot be on right edge, remapping to: "
	      << std::setprecision(dblFmt) << domainThRRad-dth << std::endl;
    srcAngleRad = domainThRRad-dth;
  }

  // find the velocity dof grid point that is closest to desired source
  auto trialDistance = std::numeric_limits<scalar_t>::max();
  for (int_t iPt=0; iPt < numGptVp; ++iPt){
    // make sure to loop over the graph and get the GID from it
    // because when we deal with sample mesh, iPt != ptGID
    // since one has more state points than actual residual points
    const auto & ptGID      = gidsList(iPt);
    const auto thisPtTh     = coords(ptGID, 0);
    const auto thisPtRadius = one/coords(ptGID, 1);

    // polar coords distance from test point to target
    auto thisD = computePolarDistance(thisPtRadius, thisPtTh, srcRadiusKm*thousand, srcAngleRad);
    if ( thisD <= trialDistance ){
      trialDistance = thisD;
      pointGid = ptGID;
    }
  }
  const auto srcTh = coords(pointGid,0);
  const auto srcRR = one/coords(pointGid,1);

  std::cout << "Source mapped to GID = "<< pointGid
	    << " with depth (km) = "
	    << std::setprecision(dblFmt) << (domainBounds[3] - srcRR)/thousand
	    << " with angle (rad) = "
	    << std::setprecision(dblFmt) << (coords(pointGid,0))
	    << " angle (deg) = "
	    << std::setprecision(dblFmt) << radToDeg(srcTh)
	    << std::endl;
}

#endif
