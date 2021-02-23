
#ifndef LOAD_BASIS_HPP_
#define LOAD_BASIS_HPP_

#include <Kokkos_Random.hpp>

namespace kokkosapp{

template <typename int_t, typename sc_t, typename basis_h_t,
	  typename parser_t, typename basis_d_t>
void loadBasis(const parser_t & parser,
	       const basis_d_t phiVp_d,
	       const basis_d_t phiSp_d)
{
  basis_h_t phiVp_h = Kokkos::create_mirror_view(phiVp_d);
  basis_h_t phiSp_h = Kokkos::create_mirror_view(phiSp_d);

  if (parser.enableRandomDummyBasis()){
    using exe_space = typename basis_h_t::execution_space;

    // dummy basis, set to random
    Kokkos::Random_XorShift64_Pool<exe_space> rand_pool(13718);
    Kokkos::fill_random(phiVp_h, rand_pool, static_cast<sc_t>(0.001));
    Kokkos::fill_random(phiSp_h, rand_pool, static_cast<sc_t>(0.001));
  }
  else{

    std::cout << "Loading basis... ";
    const auto vpUseBio  = parser.readBinaryBasis(dofId::vp);
    const auto spUseBio  = parser.readBinaryBasis(dofId::sp);
    const auto vpBasFile = parser.getBasisFileName(dofId::vp);
    const auto spBasFile = parser.getBasisFileName(dofId::sp);

    phiVp_h = readBasis<sc_t,int_t,basis_h_t>(vpBasFile, phiVp_d.extent(1), vpUseBio);
    phiSp_h = readBasis<sc_t,int_t,basis_h_t>(spBasFile, phiSp_d.extent(1), spUseBio);
    std::cout << "Done" << std::endl;
  }

  Kokkos::deep_copy(phiVp_d, phiVp_h);
  Kokkos::deep_copy(phiSp_d, phiSp_h);
}

}
#endif
