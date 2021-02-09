
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../src)

add_subdirectory(meshInfo)
add_subdirectory(parser)
add_subdirectory(seismogram)
add_subdirectory(forcing_rank1)
add_subdirectory(graphs)
add_subdirectory(coords)
add_subdirectory(jacobian_vp)
add_subdirectory(jacobian_sp)
add_subdirectory(stress_labels)

add_subdirectory(fomInnerDomain)
add_subdirectory(fomNearEarthSurface)
add_subdirectory(fomNearCmb)
add_subdirectory(fomSymmetryAxisThetaZero)
add_subdirectory(fomSymmetryAxisThetaPi)

add_subdirectory(multiDepthsForcingRank1)
add_subdirectory(multiPeriodsForcingRank1)
add_subdirectory(multiDepthsAndPeriodsForcingRank1)

