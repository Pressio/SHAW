
#ifndef ANGULAR_HELPERS_HPP_
#define ANGULAR_HELPERS_HPP_

template <typename scalar_type>
scalar_type computeCotangent(const scalar_type & theta)
{
  constexpr auto zero = constants<scalar_type>::zero();
  const auto cosine   = std::cos(theta);
  const auto sine     = std::sin(theta);

  if ( essentiallyEqual(sine, zero) )
    return zero;
  else
    return cosine/sine;
}

template <typename scalar_type>
scalar_type radToDeg(const scalar_type & valIn)
{
  constexpr auto oneEighty = static_cast<scalar_type>(180);
  constexpr auto factor = oneEighty/M_PI;
  return valIn*factor;
}

template <typename scalar_type>
scalar_type degToRad(const scalar_type & valIn)
{
  constexpr auto oneEighty = static_cast<scalar_type>(180);
  constexpr auto factor = M_PI/oneEighty;
  return valIn*factor;
}


template <typename scalar_type>
scalar_type computePolarDistance(const scalar_type & r1,
				 const scalar_type & th1,
				 const scalar_type & r2,
				 const scalar_type & th2)
{
  constexpr auto two = constants<scalar_type>::two();
  // polar coords distance from test point to target
  return std::sqrt( (r1*r1) + (r2*r2) - two*r1*r2*std::cos(th2-th1) );
}

#endif
