/*
//@HEADER
// ************************************************************************
//
// angular_helpers.hpp
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
