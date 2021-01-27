
#ifndef CHECK_CFL_HPP_
#define CHECK_CFL_HPP_

template <typename mesh_info_t, typename sc_t>
void checkCfl(const mesh_info_t & meshInfo,
	      const sc_t & dt,
	      const sc_t & maxVel) // m/s
{
  constexpr auto cfl = constants<sc_t>::cfl();
  constexpr auto two = constants<sc_t>::two();

  const auto drr     = meshInfo.getRadialSpacing(); //meters
  const auto maxArc  = meshInfo.getMaxArc(); //meters
  const auto minArc  = meshInfo.getMinArc(); //meters
  const auto hRef    = std::min(drr, std::min(maxArc, minArc));
  const auto value   = (dt*std::sqrt(two)*maxVel)/hRef;

  if (value <= cfl){
    std::cout << "CFL: OK! " << std::endl;
    std::cout << "cfl = " << value << std::endl;
  }
  else{
    std::cout << "CFL = " << value << std::endl;
    throw std::runtime_error("CFL error");
  }

}

#endif
