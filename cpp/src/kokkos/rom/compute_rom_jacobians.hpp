
#ifndef COMPUTE_ROM_JACOBIANS_HPP_
#define COMPUTE_ROM_JACOBIANS_HPP_

#include "KokkosSparse_spmv.hpp"
#include "KokkosBlas3_gemm.hpp"
#include <Kokkos_Random.hpp>

namespace kokkosapp{

template <typename parser_t, typename fom_t, typename basis_d_t, typename rom_jac_d_t>
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
    using exe_space = typename rom_jac_d_t::execution_space;

    // set zero
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

template <typename fom_t, typename basis_d_t, typename rom_jac_d_t>
void computeRomOperatorsUsingFomJacsWithoutMatProp(const fom_t & fomObj,
						   const basis_d_t & phiVp_d,
						   const basis_d_t & phiSp_d,
						   rom_jac_d_t romJvp_d,
						   rom_jac_d_t romJsp_d)
{
  using sc_t = typename fom_t::scalar_type;

  // const auto startTime = std::chrono::high_resolution_clock::now();
  // const char ct_N	= 'N';
  // const char ct_T	= 'T';
  // constexpr auto zero	= constants<sc_t>::zero();
  // constexpr auto one	= constants<sc_t>::one();

  // const auto & fomJacVp_d = fomObj.viewJacobianDevice(dofId::vp);
  // const auto & fomJacSp_d = fomObj.viewJacobianDevice(dofId::sp);
  // const auto & rhoInvVp_d = fomObj.viewInvDensityDevice(dofId::vp);
  // const auto & shearMod_d = fomObj.viewShearModulusDevice(dofId::sp);

  // rom_jac_d_t phiVpRhoInv_d("phiVpTrhoInv", nVpFom, nVp);//phiVp_d.extent(0), nVp);
  // for (std::size_t j=0; j<phiVp_d.extent(1); ++j){
  //   auto phiCol = Kokkos::subview(phiVp_d, Kokkos::ALL(), j);
  //   auto resCol = Kokkos::subview(phiVpRhoInv_d, Kokkos::ALL(), j);
  //   KokkosBlas::mult(zero, resCol, one, phiCol, rhoInvVp_d);
  // }

  // rom_jac_d_t phiSpShear_d("phiSpShear", nSpFom, nSp); //phiSp_d.extent(0), nSp);
  // for (std::size_t j=0; j<phiSp_d.extent(1); ++j){
  //   auto phiCol = Kokkos::subview(phiSp_d, Kokkos::ALL(), j);
  //   auto resCol = Kokkos::subview(phiSpShear_d, Kokkos::ALL(), j);
  //   KokkosBlas::mult(zero, resCol, one, phiCol, shearMod_d);
  // }

  // //*** do velocity ***
  // // tmp1 = fomJvp * phiSp
  // basis_d_t tmp1_d("tmp1_d", phiVp_d.extent(0), phiSp_d.extent(1));
  // KokkosSparse::spmv(&ct_N, one, fomJvp_d, phiSp_d, zero, tmp1_d);
  // // romJvp = phiVp^T * rhoVpInv * tmp1
  // KokkosBlas::gemm(&ct_T, &ct_N, one, phiVpRhoInv_d, tmp1_d, zero, romJvp_d);

  // //*** do stress ***
  // // tmp1 = fomJsp * phiVp
  // basis_d_t tmp1_d("tmp1_d", phiSp_d.extent(0), phiVp_d.extent(1));
  // KokkosSparse::spmv(&ct_N, one, fomJsp_d, phiVp_d, zero, tmp1_d);
  // // romJsp = phiSp^T * shear * tmp1
  // KokkosBlas::gemm(&ct_T, &ct_N, one, phiSpShear_d, tmp1_d, zero, romJsp_d);

  // const auto finishTime = std::chrono::high_resolution_clock::now();
  // const std::chrono::duration<double> elapsed = finishTime - startTime;
  // std::cout << "romJac time: " << std::fixed << std::setprecision(10)
  // 	    << elapsed.count() << std::endl;
}

}
#endif
