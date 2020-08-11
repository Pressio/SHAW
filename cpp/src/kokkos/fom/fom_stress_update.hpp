
#ifndef LEAP_FROG_FOM_STRESS_UPDATE_HPP_
#define LEAP_FROG_FOM_STRESS_UPDATE_HPP_

#include "KokkosSparse_spmv.hpp"
#include "KokkosBlas1_mult.hpp"

namespace kokkosapp{

template <typename sc_t, typename state_d_t, typename jac_d_t, typename shear_mod_d_t>
typename std::enable_if<
  is_accessible_on_host<state_d_t>::value and is_kokkos_1dview<state_d_t>::value
  >::type
updateStress(const bool jacobContainsMatProp,
	     const sc_t & dt,
	     state_d_t xSp_d,
	     state_d_t tmpSp,
	     const typename state_d_t::const_type xVp_d,
	     const jac_d_t jacSp_d,
	     const shear_mod_d_t shearMod_d)
{
  constexpr auto zero = constants<sc_t>::zero();
  constexpr auto one  = constants<sc_t>::one();

  if (jacobContainsMatProp){
    // xSp = xSp + dt * Jac * xVp
    KokkosSparse::spmv(KokkosSparse::NoTranspose, dt, jacSp_d, xVp_d, one, xSp_d);
  }
  else{
    // xSp = xSp + dt * (shearMod*jacSp*xVp);

    // tmpSp = Jsp * xVp
    KokkosSparse::spmv(KokkosSparse::NoTranspose, one, jacSp_d, xVp_d, zero, tmpSp);
    // xSp = xSp + dt * mu * tmpSp
    KokkosBlas::mult(one, xSp_d, dt, shearMod_d, tmpSp);
  }
}


template <typename sc_t, typename state_d_t, typename jac_d_t, typename shear_mod_d_t>
typename std::enable_if<
  is_accessible_on_host<state_d_t>::value and is_kokkos_2dview<state_d_t>::value
  >::type
updateStress(const bool jacobContainsMatProp,
	     const sc_t & dt,
	     state_d_t xSp_d,
	     state_d_t tmpSp,
	     const typename state_d_t::const_type xVp_d,
	     const jac_d_t jacSp_d,
	     const shear_mod_d_t shearMod_d)
{
  constexpr auto one  = constants<sc_t>::one();

  // xSp = xSp + dt * Jac * xVp
  KokkosSparse::spmv(KokkosSparse::NoTranspose, dt, jacSp_d, xVp_d, one, xSp_d);

}

}//end namespace kokkosapp
#endif
