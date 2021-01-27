
#ifndef COMPUTE_NUM_VP_POINTS_ON_SURFACE_HPP_
#define COMPUTE_NUM_VP_POINTS_ON_SURFACE_HPP_

#include <iostream>

template <typename sc_t, typename int_t, typename gids_t, typename coords_t>
void countVpPointsOnEarthSurface(const int_t numGptVp,
				 const gids_t & gidsVp,
				 const coords_t & coordsVp,
				 const sc_t & earthSurfRadius, //[m]
				 int_t & result)
{
  constexpr auto one	= constants<sc_t>::one();

  result = constants<int_t>::zero();
  for (int_t iPt=0; iPt < numGptVp; ++iPt){
    const auto & ptGID	= gidsVp(iPt);
    const auto rr	= one/coordsVp(ptGID, 1);
    if ( essentiallyEqual(rr, earthSurfRadius) )
      result++;
  }
  std::cout << "NumVpPts on surface  = " << result << std::endl;
}

#endif
