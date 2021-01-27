
#ifndef RANK_ONE_FORCING_HPP_
#define RANK_ONE_FORCING_HPP_

#include "map_point_source_to_velocity_grid_point.hpp"
#include "KokkosBlas1_fill.hpp"
#include "KokkosBlas1_scal.hpp"

template <typename sc_t, typename state_d_t, typename int_t>
class RankOneForcing
{
  static_assert(is_vector_kokkos< state_d_t >::value,
		"Rank-1 forcing must use a rank-1 kokkos view");

  using state_h_t = typename state_d_t::host_mirror_type;

  // f_h_ contains the full time series of the signal
  state_h_t f_h_;

  // f_d_ contains the forcing vector over the mesh
  state_d_t f_d_;
  // // indicator_d_ is all zeros except for the entry where forcing acts
  // // this is set upon construction, and used to evaluate the forcing
  // // each time the evaluation is requested, see evaluate method below
  // state_d_t indicator_d_;

  // The signal is always associated with a velocity point, not a stress point.
  // myVpGid_ identifies which velocity point the signal is located at
  int_t myVpGid_ = -1;

  sc_t maxFreq_ = {};
  const sc_t dt_ = {};
  const int_t NSteps_ = {};

public:
  template <typename signal_t, typename parser_t, typename mesh_info_t, typename app_t>
  RankOneForcing(const signal_t & signalObj,
		 const parser_t & parser,
		 const mesh_info_t & meshInfo,
		 const app_t & appObj,
		 const sc_t depthKm,
		 const sc_t angleDeg)
    : f_h_("Fh", parser.getNumSteps()),
      f_d_("Fd", meshInfo.getNumResidualVpPts()),
      dt_(parser.getTimeStepSize()),
      NSteps_(parser.getNumSteps()),
      /*indicator_d_("Indic_d", meshInfo.getNumResidualVpPts()),*/
      maxFreq_(signalObj.getFrequency())
  {
    const auto gidsVp = appObj.viewGidListHost(dofId::vp);
    const auto coords = appObj.viewCoordsHost(dofId::vp);

    // find the vpGid identifying the grid point where the source is mapped to
    const sc_t myRadiusKm = constants<sc_t>::earthSurfaceRadiusKm() - depthKm;
    mapPointSourceToGridPoint(angleDeg, myRadiusKm, depthKm,
			      meshInfo.viewDomainBounds(),
			      meshInfo.getNumResidualVpPts(), gidsVp, coords,
			      meshInfo.getAngularSpacing(), myVpGid_);

    KokkosBlas::fill(f_h_, constants<sc_t>::zero());
    KokkosBlas::fill(f_d_, constants<sc_t>::zero());

    // // set the target entry in the indicator host vector equal to 1
    // state_h_t indicator_h = Kokkos::create_mirror_view(indicator_d_);
    // KokkosBlas::fill(indicator_h, constants<sc_t>::zero());
    // indicator_h(myVpGid_) = constants<sc_t>::one();
    // Kokkos::deep_copy(indicator_d_, indicator_h);

    // store the full time series of the signal into the host array
    storeSignalTimeSeries(signalObj);
  }

  template <typename parser_t, typename mesh_info_t, typename app_t>
  RankOneForcing(const parser_t    & parser,
		 const mesh_info_t & meshInfo,
		 const app_t	   & appObj)
    : RankOneForcing(Signal<sc_t>(parser.getSourceSignalKind(),
				  parser.getSourceProperty("delay"),
				  parser.getSourceProperty("period")),
		     parser, meshInfo, appObj,
		     parser.getSourceProperty("depth"),
		     parser.getSourceProperty("angle"))
  {}

  template <typename signal_t>
  void replaceSignal(const signal_t & newSignal)
  {
    maxFreq_ = newSignal.getFrequency();
    KokkosBlas::fill(f_h_, constants<sc_t>::zero());
    storeSignalTimeSeries(newSignal);
  }

private:
  template <typename signal_t>
  void storeSignalTimeSeries(const signal_t signal){
    // store the full time series of the signal into the host array
    sc_t time = constants<sc_t>::zero();
    for (int_t iStep = 1; iStep<=NSteps_; ++iStep){
      signal(time, f_h_(iStep-1));
      time = iStep * dt_;
    }
  }

public:
  // for single source, max frquency is the frequency of the source signal
  sc_t getMaxFreq() const{ return maxFreq_; }

  int_t getVpGid() const{ return myVpGid_; }

  sc_t getForcingValueAtStep(std::size_t step) const{
    return f_h_(step-1);
  }

  state_d_t viewForcingDevice() const{ return f_d_; }

  void evaluate(sc_t time, std::size_t step){
    KokkosBlas::fill(f_d_, constants<sc_t>::zero());
    const auto src = Kokkos::subview(f_h_, step-1);
    const auto des = Kokkos::subview(f_d_, myVpGid_);
    Kokkos::deep_copy(des, src);
  }

  void complexityOfEvaluateMethod(double & memCostMB, double flopsCost) const{
    // no operation is done during evaluate, just copying, see above
    const double memMBCostFill = 1.*( f_d_.extent(0)*sizeof(sc_t) )/1024./1024.;
    // for copy we have one read + one write
    const double memMBCostCopy = 1.*( 2*sizeof(sc_t) )/1024./1024.;

    memCostMB = memMBCostFill + memMBCostCopy;
    flopsCost = 0.;
  }

  // // // *** for non-host ****
  // // template <typename _state_d_t = state_d_t>
  // // typename std::enable_if< !is_accessible_on_host<_state_d_t>::value >::type
  // // evaluate(sc_t time, std::size_t step){
  // //   KokkosBlas::scal(f_d_, f_h_(step-1), indicator_d_);
  // // }

  // template <typename _state_d_t = state_d_t>
  // typename std::enable_if< !is_accessible_on_host<_state_d_t>::value >::type
  // complexity(double & memCostMB, double flopsCost) const{
  //   using ord_t = typename state_d_t::traits::size_type;
  //   using comp_t = Complexity<sc_t, ord_t>;
  //   comp_t::scal(f_d_.extent(0), memCostMB, flopsCost);
  // }
};

#endif
