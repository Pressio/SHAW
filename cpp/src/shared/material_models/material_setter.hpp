
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
