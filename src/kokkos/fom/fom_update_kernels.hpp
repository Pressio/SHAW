
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

  // maybe we should do the following on host directly since
  // for a single forcing, if pointwise, we only change a single element
  auto f_d = fObj.viewForcingDevice();
  KokkosBlas::mult(one, xVp_d, dt, rhoInvVp_d, f_d);
  Kokkos::fence();
}

// rank-1 specialize
template <typename sc_t, typename state_d_t, typename jac_d_t>
typename std::enable_if<is_kokkos_1dview<state_d_t>::value>::type
updateStress(const sc_t & dt,
	     state_d_t xSp_d,
	     const typename state_d_t::const_type xVp_d,
	     jac_d_t jacSp_d)
{
  // xSp = xSp + dt * Jac * xVp

  constexpr auto one  = constants<sc_t>::one();
  KokkosSparse::spmv(KokkosSparse::NoTranspose, dt, jacSp_d, xVp_d, one, xSp_d);
  Kokkos::fence();
}

template <class sc_t, class state_t, class gids_t, class f_t, class rho_inv_t>
struct AddForcingRank2
{
  sc_t dt_;
  state_t x_;
  gids_t gids_;
  f_t f_;
  rho_inv_t rhoInv_;

  AddForcingRank2(const sc_t & dt,
		  state_t x,
		  gids_t gids,
		  f_t f,
		  rho_inv_t rhoInv)
    : dt_(dt), x_(x), gids_(gids), f_(f), rhoInv_(rhoInv){}

  KOKKOS_INLINE_FUNCTION
  void operator() (std::size_t i) const
  {
    const auto gidValue = gids_(i);
    const auto & rhoInv = rhoInv_(gidValue);
    x_(gidValue, i) += rhoInv*f_(i)*dt_;
  }
};


// rank-2 specialize
template <
  typename sc_t,
  typename state_d_t,
  typename jac_d_t,
  typename rho_inv_d_t,
  typename forcing_t
  >
typename std::enable_if<is_kokkos_2dview<state_d_t>::value>::type
updateVelocity(const sc_t & dt,
	       state_d_t xVp_d,
	       typename state_d_t::const_type xSp_d,
	       const jac_d_t jacVp_d,
	       const rho_inv_d_t rhoInvVp_d,
	       forcing_t & fObj)
{
  /*
   *	A1	xVp = xVp + dt * Jvp * xSp
   *	A2	xVp = xVp + dt * rhoInvVp * f
   */

  constexpr auto one  = constants<sc_t>::one();
  KokkosSparse::spmv(KokkosSparse::NoTranspose, dt, jacVp_d, xSp_d, one, xVp_d);
  // auto f_d = fObj.viewForcingDevice();
  // KokkosBlas::mult(one, xVp_d, dt, rhoInvVp_d, f_d);

  // maybe we should do the following on host directly since
  // it might be too small for device
  auto vpGids_d = fObj.getVpGidsDevice();
  auto f_d = fObj.viewForcingDevice();
  using gids_t = decltype(vpGids_d);
  using f_d_t  = decltype(f_d);
  using functor_t = AddForcingRank2<sc_t, state_d_t, gids_t, f_d_t, rho_inv_d_t>;
  functor_t fnc(dt, xVp_d, vpGids_d, f_d, rhoInvVp_d);
  Kokkos::parallel_for(vpGids_d.extent(0), fnc);

  Kokkos::fence();


  //constexpr auto one  = constants<sc_t>::one();
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
  //throw std::runtime_error("FOM velo update for rank-2 not implemented yet");
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
