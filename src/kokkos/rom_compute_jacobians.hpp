/*
//@HEADER
// ************************************************************************
//
// compute_rom_jacobians.hpp
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

#ifndef COMPUTE_ROM_JACOBIANS_HPP_
#define COMPUTE_ROM_JACOBIANS_HPP_

#include "KokkosSparse_spmv.hpp"
#include "KokkosBlas3_gemm.hpp"
#include <Kokkos_Random.hpp>

namespace kokkosapp{

template <class parser_t, class fom_t, class basis_d_t, class rom_jac_d_t>
void computeRomOperatorsUsingFomJacsWithMatProp(const parser_t & parser,
						const fom_t & fomObj,
						const basis_d_t phiVp_d,
						const basis_d_t phiSp_d,
						rom_jac_d_t romJvp_d,
						rom_jac_d_t romJsp_d)
{
  using sc_t = typename fom_t::scalar_type;

  const auto startTime = std::chrono::high_resolution_clock::now();
  if (parser.disableCompRomJacobians())
  {
    constexpr auto zero	= constants<sc_t>::zero();
    KokkosBlas::fill(romJvp_d, zero);
    KokkosBlas::fill(romJsp_d, zero);
    std::cout << "The computation of the rom jacobians is disabled" << std::endl;
  }
  else{

    const char ct_N	= 'N';
    const char ct_T	= 'T';
    constexpr auto zero	= constants<sc_t>::zero();
    constexpr auto one	= constants<sc_t>::one();

    const auto fomJvp_d = fomObj.viewJacobianDevice(dofId::vp);
    const auto fomJsp_d = fomObj.viewJacobianDevice(dofId::sp);

    //*** do velocity ***
    {
      // tmp1 = fomJvp * phiSp
      basis_d_t tmp1_d("tmp1_d", phiVp_d.extent(0), phiSp_d.extent(1));
      KokkosSparse::spmv(&ct_N, one, fomJvp_d, phiSp_d, zero, tmp1_d);
      // romJvp = phiVp^T * tmp1
      KokkosBlas::gemm(&ct_T, &ct_N, one, phiVp_d, tmp1_d, zero, romJvp_d);
    }

    //*** do stress ***
    {
      // tmp1 = fomJsp * phiVp
      basis_d_t tmp1_d("tmp1_d", phiSp_d.extent(0), phiVp_d.extent(1));
      KokkosSparse::spmv(&ct_N, one, fomJsp_d, phiVp_d, zero, tmp1_d);
      // romJsp = phiSp^T * tmp1
      KokkosBlas::gemm(&ct_T, &ct_N, one, phiSp_d, tmp1_d, zero, romJsp_d);
    }
  }

  const auto finishTime = std::chrono::high_resolution_clock::now();
  const std::chrono::duration<double> elapsed = finishTime - startTime;
  std::cout << "romJac time: " << std::fixed << std::setprecision(10)
	    << elapsed.count() << std::endl;
}

}
#endif


// template <typename fom_t, typename basis_d_t, typename rom_jac_d_t>
// void computeRomOperatorsUsingFomJacsWithoutMatProp(const fom_t & fomObj,
// 						   const basis_d_t & phiVp_d,
// 						   const basis_d_t & phiSp_d,
// 						   rom_jac_d_t romJvp_d,
// 						   rom_jac_d_t romJsp_d)
// {
//   using sc_t = typename fom_t::scalar_type;

//   // const auto startTime = std::chrono::high_resolution_clock::now();
//   // const char ct_N	= 'N';
//   // const char ct_T	= 'T';
//   // constexpr auto zero	= constants<sc_t>::zero();
//   // constexpr auto one	= constants<sc_t>::one();

//   // const auto & fomJacVp_d = fomObj.viewJacobianDevice(dofId::vp);
//   // const auto & fomJacSp_d = fomObj.viewJacobianDevice(dofId::sp);
//   // const auto & rhoInvVp_d = fomObj.viewInvDensityDevice(dofId::vp);
//   // const auto & shearMod_d = fomObj.viewShearModulusDevice(dofId::sp);

//   // rom_jac_d_t phiVpRhoInv_d("phiVpTrhoInv", nVpFom, nVp);//phiVp_d.extent(0), nVp);
//   // for (std::size_t j=0; j<phiVp_d.extent(1); ++j){
//   //   auto phiCol = Kokkos::subview(phiVp_d, Kokkos::ALL(), j);
//   //   auto resCol = Kokkos::subview(phiVpRhoInv_d, Kokkos::ALL(), j);
//   //   KokkosBlas::mult(zero, resCol, one, phiCol, rhoInvVp_d);
//   // }

//   // rom_jac_d_t phiSpShear_d("phiSpShear", nSpFom, nSp); //phiSp_d.extent(0), nSp);
//   // for (std::size_t j=0; j<phiSp_d.extent(1); ++j){
//   //   auto phiCol = Kokkos::subview(phiSp_d, Kokkos::ALL(), j);
//   //   auto resCol = Kokkos::subview(phiSpShear_d, Kokkos::ALL(), j);
//   //   KokkosBlas::mult(zero, resCol, one, phiCol, shearMod_d);
//   // }

//   // //*** do velocity ***
//   // // tmp1 = fomJvp * phiSp
//   // basis_d_t tmp1_d("tmp1_d", phiVp_d.extent(0), phiSp_d.extent(1));
//   // KokkosSparse::spmv(&ct_N, one, fomJvp_d, phiSp_d, zero, tmp1_d);
//   // // romJvp = phiVp^T * rhoVpInv * tmp1
//   // KokkosBlas::gemm(&ct_T, &ct_N, one, phiVpRhoInv_d, tmp1_d, zero, romJvp_d);

//   // //*** do stress ***
//   // // tmp1 = fomJsp * phiVp
//   // basis_d_t tmp1_d("tmp1_d", phiSp_d.extent(0), phiVp_d.extent(1));
//   // KokkosSparse::spmv(&ct_N, one, fomJsp_d, phiVp_d, zero, tmp1_d);
//   // // romJsp = phiSp^T * shear * tmp1
//   // KokkosBlas::gemm(&ct_T, &ct_N, one, phiSpShear_d, tmp1_d, zero, romJsp_d);

//   // const auto finishTime = std::chrono::high_resolution_clock::now();
//   // const std::chrono::duration<double> elapsed = finishTime - startTime;
//   // std::cout << "romJac time: " << std::fixed << std::setprecision(10)
//   // 	    << elapsed.count() << std::endl;
// }
