
#ifndef LEAP_FROG_RUN_FOM_HPP_
#define LEAP_FROG_RUN_FOM_HPP_

#include "fom_velocity_update.hpp"
#include "fom_stress_update.hpp"
#include "fom_complexities.hpp"
#include <numeric>

namespace kokkosapp{

template<typename state_d_t>
typename std::enable_if< is_kokkos_1dview<state_d_t>::value, state_d_t >::type
createAuxiliaryStateIfNeeded(state_d_t x, bool jacobContainsMatProp){
  state_d_t auxX("auxX", 1);
  if (!jacobContainsMatProp) Kokkos::resize(auxX, x.extent(0));
  return auxX;
}

template<typename state_d_t>
typename std::enable_if< is_kokkos_2dview<state_d_t>::value, state_d_t >::type
createAuxiliaryStateIfNeeded(state_d_t x, bool jacobContainsMatProp){
  state_d_t auxX("auxX", 1, 1);
  if (!jacobContainsMatProp) Kokkos::resize(auxX, x.extent(0), x.extent(1));
  return auxX;
}


template <
  typename step_t, typename sc_t, typename app_t, typename forcing_t,
  typename observer_t, typename seismo_t, typename state_d_t
  >
void runFom(bool jacobContainsMatProp,
	    bool exploitForcingSparsity,
	    const step_t & numSteps,
	    const sc_t dt,
	    const app_t & fomObj,
	    forcing_t & forcingObj,
	    observer_t & observerObj,
	    seismo_t & seismoObj,
	    state_d_t xVp_d,
	    state_d_t xSp_d)
{
  KokkosBlas::fill(xVp_d, constants<sc_t>::zero());
  KokkosBlas::fill(xSp_d, constants<sc_t>::zero());

  const auto jacVp_d     = fomObj.viewJacobianDevice(dofId::vp);
  const auto jacSp_d     = fomObj.viewJacobianDevice(dofId::sp);
  const auto rhoInvVp_d  = fomObj.viewInvDensityDevice(dofId::vp);
  const auto shearMod_d  = fomObj.viewShearModulusDevice(dofId::sp);

  // xSpAux_d is only needed when the material properties are factored out the Jacobians
  // so let it initialized to 1 entry, and resize appropriately if needed
  state_d_t xSpAux_d = createAuxiliaryStateIfNeeded(xSp_d, jacobContainsMatProp);

  // to collec timings
  Kokkos::Timer timer;
  double dataCollectionTime = {};
  std::array<double, 3> perfTimes = {1e32,0.,0.}; //min, max, total

  //****** LOOP ******//
  const auto startTime  = std::chrono::high_resolution_clock::now();
  sc_t timeVp = constants<sc_t>::zero();
  for (auto iStep = 1; iStep<=numSteps; ++iStep)
  {
    if (iStep % 2000 == 0) std::cout << "FOM-KOKKOS step = " << iStep << std::endl;

    // compute forcing for current time
    timer.reset();
    forcingObj.evaluate(timeVp, iStep);
    const double ct1 = timer.seconds();

    // ----------------
    // 1. do velocity
    timer.reset();
    updateVelocity(jacobContainsMatProp, exploitForcingSparsity,
		   dt, iStep, timeVp, xVp_d, xSp_d,
		   jacVp_d, rhoInvVp_d, forcingObj);
    const double ct2 = timer.seconds();
    timer.reset();
    observerObj.observe(dofId::vp, iStep, xVp_d);
    seismoObj.storeVelocitySignalAtReceivers(iStep, xVp_d);
    dataCollectionTime += timer.seconds();

    // update time
    timeVp = iStep*dt;

    // ----------------
    // 2. do stress
    timer.reset();
    updateStress(jacobContainsMatProp, dt, xSp_d, xSpAux_d, xVp_d, jacSp_d, shearMod_d);
    const double ct3 = timer.seconds();
    timer.reset();
    observerObj.observe(dofId::sp, iStep, xSp_d);
    dataCollectionTime += timer.seconds();

    // ----------------
    // 3. timing vars
    const double time = ct1+ct2+ct3;
    perfTimes[0] = std::min(perfTimes[0], time);
    perfTimes[1] = std::max(perfTimes[1], time);
    perfTimes[2] += time;
  }

  const auto finishTime = std::chrono::high_resolution_clock::now();
  const std::chrono::duration<double> elapsed = finishTime - startTime;
  std::cout << "\nloopTime = " << std::fixed << std::setprecision(10) << elapsed.count();
  std::cout << "\ndataCollectionTime = " << std::fixed << std::setprecision(10)
	    << dataCollectionTime << std::endl;

  // compute complexity and print
  double memCostMB, flopsCost = 0.;
  complexityFom<sc_t>(jacobContainsMatProp, exploitForcingSparsity,
		      xVp_d, xSp_d, fomObj, forcingObj, memCostMB, flopsCost);
  printPerf(numSteps, perfTimes, memCostMB, flopsCost);
}

}//end namespace kokkosapp
#endif
