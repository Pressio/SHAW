
#ifndef CHECK_DISPERSION_CRITERION_HPP_
#define CHECK_DISPERSION_CRITERION_HPP_

template <typename mesh_info_t, typename sc_t>
void checkDispersionCriterion(const mesh_info_t & meshInfo,
			      const sc_t & maxfreq,
			      const sc_t & minVel) // m/s
{
  const auto drr	 = meshInfo.getRadialSpacing(); //meters
  const auto minArc	 = meshInfo.getMinArc();	//meters
  const auto maxArc	 = meshInfo.getMaxArc();	//meters
  const auto ratio	 = minVel/((sc_t) Nlambda * maxfreq);

  const auto f0 = minVel / ((sc_t) Nlambda * std::max(drr, maxArc));
  std::cout << "targetPeriod = " << 1./maxfreq << std::endl;
  std::cout << "centerFreq = " << f0 << std::endl;
  std::cout << "minPeriod  = " << 1./f0 << std::endl;

  if (drr <= ratio){
    std::cout << "Numerical dispersion criterion along r: OK! " << std::endl;
    std::cout << "drr = " << drr << ", ratio = " << ratio << std::endl;
  }
  else{
    std::cout << "drr = " << drr << ", ratio = " << ratio << std::endl;
    throw std::runtime_error("Numerical dispersion criterion violated by radial spacing");
  }

  if (maxArc <= ratio){
    std::cout << "Numerical dispersion criterion along theta: OK! " << std::endl;
    std::cout << "maxArc = " << maxArc << ", ratio = " << ratio << std::endl;
  }
  else{
    std::cout << "maxArc = " << maxArc << ", ratio = " << ratio << std::endl;
    throw std::runtime_error("Numerical dispersion criterion violated along theta");
  }
}

#endif
