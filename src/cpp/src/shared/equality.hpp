
#ifndef EQUALITY_HPP_
#define EQUALITY_HPP_

#include <limits>

template <typename sc_t>
bool essentiallyEqual(sc_t a, sc_t b){
  constexpr auto eps = std::numeric_limits<sc_t>::epsilon();
  return std::abs(a - b) <= ( (std::abs(a) > std::abs(b) ? std::abs(b) : std::abs(a)) * eps);
}

#endif
