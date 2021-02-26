/*
//@HEADER
// ************************************************************************
//
// check_cfl.hpp
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
