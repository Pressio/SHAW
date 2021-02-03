
#ifndef LEAP_FROG_FOM_VELOCITY_UPDATE_HPP_
#define LEAP_FROG_FOM_VELOCITY_UPDATE_HPP_

#include "KokkosSparse_spmv.hpp"
#include "KokkosBlas1_mult.hpp"
#include "KokkosBlas1_axpby.hpp"
#include <numeric>

namespace kokkosapp{

// rank-1 specialize
template <
  typename sc_t,
  typename state_d_t,
  typename jac_d_t,
  typename rho_inv_d_t,
  typename forcing_t
  >
typename std::enable_if<is_kokkos_1dview<state_d_t>::value>::type
updateVelocity(const sc_t & dt,
	       state_d_t xVp_d,
	       typename state_d_t::const_type xSp_d,
	       const jac_d_t jacVp_d,
	       const rho_inv_d_t rhoInvVp_d,
	       forcing_t & fObj)
{
  /* compute the velocity update:
   *	xVp = xVp + dt*jacVp*xSp + dt*rhoInvVp*f;
   *
   * which we split into two steps:
   *	A1	xVp = xVp + dt * Jvp * xSp
   *	A2	xVp = xVp + dt * rhoInvVp * f
   */

  constexpr auto one  = constants<sc_t>::one();

  KokkosSparse::spmv(KokkosSparse::NoTranspose, dt, jacVp_d, xSp_d, one, xVp_d);

  auto f_d = fObj.viewForcingDevice();
  KokkosBlas::mult(one, xVp_d, dt, rhoInvVp_d, f_d);
}

// rank-1 specialize
template <typename sc_t, typename state_d_t, typename jac_d_t>
typename std::enable_if<is_kokkos_1dview<state_d_t>::value>::type
updateStress(const sc_t & dt,
	     state_d_t xSp_d,
	     const typename state_d_t::const_type xVp_d,
	     const jac_d_t jacSp_d)
{
  // xSp = xSp + dt * Jac * xVp

  constexpr auto zero = constants<sc_t>::zero();
  constexpr auto one  = constants<sc_t>::one();
  KokkosSparse::spmv(KokkosSparse::NoTranspose, dt, jacSp_d, xVp_d, one, xSp_d);
}

// rank-2 specialize
template <
  typename sc_t,
  typename state_d_t,
  typename jac_d_t,
  typename rho_inv_d_t,
  typename forcing_t
  >
typename std::enable_if<
  is_accessible_on_host<state_d_t>::value and is_kokkos_2dview<state_d_t>::value
  >::type
updateVelocity(const sc_t & dt,
	       state_d_t xVp_d,
	       typename state_d_t::const_type xSp_d,
	       const jac_d_t jacVp_d,
	       const rho_inv_d_t rhoInvVp_d,
	       forcing_t & fObj)
{
  constexpr auto one  = constants<sc_t>::one();

  // if (exploitForcingSparsity){
  //   // xVp = xVp + dt * Jvp * xSp
  //   KokkosSparse::spmv(KokkosSparse::NoTranspose, dt, jacVp_d, xSp_d, one, xVp_d);

  //   const auto vpGid   = fObj.getVpGid();
  //   auto f_v = fObj.getForcingAtStep(step);
  //   for (std::size_t j=0; j<xVp_d.extent(1); ++j)
  //     xVp_d(vpGid, j) += rhoInvVp_d(vpGid)*f_v(j)*dt;

  //   //auto alpha = rhoInvVp_d(vpGid)*dt;
  //   // auto x_v = Kokkos::subview(xVp_d, vpGid, Kokkos::ALL());
  //   // KokkosBlas::axpy(alpha, f_v, x_v);
  // }
  // else{
  throw std::runtime_error("FOM velo update for rank-2 not implemented yet");
  //}
}

// rank-2 specialize
template <typename sc_t, typename state_d_t, typename jac_d_t>
typename std::enable_if<is_kokkos_2dview<state_d_t>::value>::type
updateStress(const sc_t & dt,
	     state_d_t xSp_d,
	     const typename state_d_t::const_type xVp_d,
	     const jac_d_t jacSp_d)
{
  // xSp = xSp + dt * Jac * xVp

  constexpr auto one  = constants<sc_t>::one();
  KokkosSparse::spmv(KokkosSparse::NoTranspose, dt, jacSp_d, xVp_d, one, xSp_d);
}

}//end namespace kokkosapp
#endif
