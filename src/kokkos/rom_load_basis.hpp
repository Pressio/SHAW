/*
//@HEADER
// ************************************************************************
//
// load_basis.hpp
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

#ifndef LOAD_BASIS_HPP_
#define LOAD_BASIS_HPP_

#include <Kokkos_Random.hpp>

namespace kokkosapp{

template <
  class int_t,
  class sc_t,
  class parser_t,
  class basis_d_t
  >
void fillBasis(const parser_t & parser,
	       const basis_d_t phiVp_d,
	       const basis_d_t phiSp_d)
{
  auto phiVp_h = Kokkos::create_mirror_view(phiVp_d);
  auto phiSp_h = Kokkos::create_mirror_view(phiSp_d);

  if (parser.enableRandomDummyBasis()){
    using exe_space = typename decltype(phiVp_h)::execution_space;

    // dummy basis, set to random
    Kokkos::Random_XorShift64_Pool<exe_space> rand_pool(13718);
    Kokkos::fill_random(phiVp_h, rand_pool, static_cast<sc_t>(0.001));
    Kokkos::fill_random(phiSp_h, rand_pool, static_cast<sc_t>(0.001));
  }
  else{
    std::cout << "Loading basis... \n";

    const auto vpUseBio  = parser.readBinaryBasis(dofId::vp);
    const auto spUseBio  = parser.readBinaryBasis(dofId::sp);
    readBasis(parser.getBasisFileName(dofId::vp), phiVp_h, vpUseBio);
    readBasis(parser.getBasisFileName(dofId::sp), phiSp_h, spUseBio);

    double min ={};
    double max = {};
    for (std::size_t i=0; i<phiVp_h.extent(0); ++i){
      for (std::size_t j=0; j<phiVp_h.extent(1); ++j){
	min = std::min(min, phiVp_h(i,j));
	max = std::max(max, phiVp_h(i,j));
      }
    }
    std::cout << " MIN/MAX = " << min << " " << max << std::endl;


    std::cout << "Done" << std::endl;
  }

  Kokkos::deep_copy(phiVp_d, phiVp_h);
  Kokkos::deep_copy(phiSp_d, phiSp_h);
}

}
#endif
