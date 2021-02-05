
#ifndef UTILS_CONSTANTS_HPP_
#define UTILS_CONSTANTS_HPP_

#include <limits>
#include <chrono>
#include <string>
#include <cmath>

constexpr auto dblFmt = std::numeric_limits<double>::max_digits10;

// number of grid points per wavelength
static constexpr int Nlambda = 10;

template<typename scalar_t = double>
struct constants
{
  static constexpr scalar_t cfl(){ return static_cast<scalar_t>(0.28); }

  static constexpr scalar_t negOne(){ return static_cast<scalar_t>(-1); }
  static constexpr scalar_t zero(){ return static_cast<scalar_t>(0); }
  static constexpr scalar_t one(){ return static_cast<scalar_t>(1); }
  static constexpr scalar_t two(){ return static_cast<scalar_t>(2); }
  static constexpr scalar_t three(){ return static_cast<scalar_t>(3); }
  static constexpr scalar_t four(){return static_cast<scalar_t>(4);}
  static constexpr scalar_t thousand(){return static_cast<scalar_t>(1000);}
};

#endif
