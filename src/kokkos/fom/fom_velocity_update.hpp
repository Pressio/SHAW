
#ifndef LEAP_FROG_FOM_VELOCITY_UPDATE_HPP_
#define LEAP_FROG_FOM_VELOCITY_UPDATE_HPP_

#include "KokkosSparse_spmv.hpp"
#include "KokkosBlas1_mult.hpp"
#include "KokkosBlas1_axpby.hpp"
#include <numeric>

namespace kokkosapp{

template <
  typename sc_t,
  typename step_t,
  typename state_d_t,
  typename jac_d_t,
  typename rho_inv_d_t,
  typename forcing_t
  >
typename std::enable_if<
  is_accessible_on_host<state_d_t>::value and is_kokkos_1dview<state_d_t>::value
  >::type
updateVelocity(const bool includeMatPropInJac,
	       const bool exploitForcingSparsity,
	       const sc_t & dt,
	       const step_t & step,
	       const sc_t & time,
	       state_d_t xVp_d,
	       typename state_d_t::const_type xSp_d,
	       const jac_d_t jacVp_d,
	       const rho_inv_d_t rhoInvVp_d,
	       forcing_t & fObj)
{
  constexpr auto one  = constants<sc_t>::one();

  /* to compute the velocity update, we differentiate:
   * Case A: the material prop are included in the Jacobians, then we do:
   *	A1	xVp = xVp + dt * Jvp * xSp
   *	A2	xVp = xVp + dt * rhoInvVp * f
   *		result: xVp = xVp + dt*jacVp*xSp + dt*rhoInvVp*f;
   *
   * Case B: the material prop are NOT included in the Jacobians, then we do:
   *	B1	f = f + Jvp * xSp, which is fine because f is updated at every iteration
   *	B2	xVp = xVp + dt * rhoInv * f
   *		result: xVp = xVp + dt*rhoInvVp*jacVp*xSp + dt*rhoInvVp*f;
   */

  if (includeMatPropInJac){
    // A1
    KokkosSparse::spmv(KokkosSparse::NoTranspose, dt, jacVp_d, xSp_d, one, xVp_d);
    // A2
    if (exploitForcingSparsity){
      const auto vpGid   = fObj.getVpGid();
      const auto fValue  = fObj.getForcingValueAtStep(step);
      xVp_d(vpGid) += dt*rhoInvVp_d(vpGid)*fValue;
    }else{
      auto f_d = fObj.viewForcingDevice();
      KokkosBlas::mult(one, xVp_d, dt, rhoInvVp_d, f_d);
    }
  }
  else{
    //B1
    auto f_d = fObj.viewForcingDevice();
    KokkosSparse::spmv(KokkosSparse::NoTranspose, one, jacVp_d, xSp_d, one, f_d);
    //B2
    KokkosBlas::mult(one, xVp_d, dt, rhoInvVp_d, f_d);
  }
}

template <
  typename sc_t,
  typename step_t,
  typename state_d_t,
  typename jac_d_t,
  typename rho_inv_d_t,
  typename forcing_t
  >
typename std::enable_if<
  is_accessible_on_host<state_d_t>::value and is_kokkos_2dview<state_d_t>::value
  >::type
updateVelocity(const bool includeMatPropInJac,
	       const bool exploitForcingSparsity,
	       const sc_t & dt,
	       const step_t & step,
	       const sc_t & time,
	       state_d_t xVp_d,
	       typename state_d_t::const_type xSp_d,
	       const jac_d_t jacVp_d,
	       const rho_inv_d_t rhoInvVp_d,
	       forcing_t & fObj)
{
  constexpr auto one  = constants<sc_t>::one();

  if (exploitForcingSparsity){
    // xVp = xVp + dt * Jvp * xSp
    KokkosSparse::spmv(KokkosSparse::NoTranspose, dt, jacVp_d, xSp_d, one, xVp_d);

    const auto vpGid   = fObj.getVpGid();
    auto f_v = fObj.getForcingAtStep(step);
    for (std::size_t j=0; j<xVp_d.extent(1); ++j)
      xVp_d(vpGid, j) += rhoInvVp_d(vpGid)*f_v(j)*dt;

    //auto alpha = rhoInvVp_d(vpGid)*dt;
    // auto x_v = Kokkos::subview(xVp_d, vpGid, Kokkos::ALL());
    // KokkosBlas::axpy(alpha, f_v, x_v);
  }
  else{
    throw std::runtime_error("Calling impl for FOM with rank-2 non-sparse forcing not impl yet");
  }
}

}//end namespace kokkosapp
#endif










// template <typename sc_t, typename state_d_t>
// struct Foo2
// {
//   std::size_t vpGid_;
//   sc_t T_;
//   sc_t dt_;
//   state_d_t xVp_d_;
//   state_d_t rhoInvVp_d_;

//   Foo2(sc_t T, sc_t dt, state_d_t xVp_d, state_d_t rhoInvVp, std::size_t vpGid)
//     : T_(T), dt_(dt), xVp_d_(xVp_d), rhoInvVp_d_(rhoInvVp), vpGid_(vpGid)
//   {}

//   KOKKOS_INLINE_FUNCTION
//   void operator() (const std::size_t & i) const
//   {
//     const sc_t freq = 1./50.;
//     const sc_t delayTime = 120.;

//     const auto tDiff   = (T_-delayTime);
//     const auto tDiffSq = tDiff*tDiff;
//     const auto expTerm = std::exp( -freq*tDiffSq );
//     const auto result = (i == vpGid_) ? -2.*tDiff*freq*expTerm : 0.0;
//     xVp_d_(i) += dt_*rhoInvVp_d_(i) * result;
//   }
// };
