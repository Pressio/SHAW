
#ifndef RANK_TWO_FORCING_HPP_
#define RANK_TWO_FORCING_HPP_

#include "KokkosBlas1_fill.hpp"
#include "KokkosBlas1_scal.hpp"

template <typename sc_t, typename signal_instances_h_type>
class RankTwoForcing
{
public:
  using f_d_t = Kokkos::View<sc_t*, Kokkos::DefaultExecutionSpace>;
  using f_h_t = typename f_d_t::host_mirror_type;
  using vp_gids_d_t = Kokkos::View<std::size_t*, Kokkos::DefaultExecutionSpace>;
  using vp_gids_h_t = typename vp_gids_d_t::host_mirror_type;

private:
  // list of signal objects
  signal_instances_h_type signals_;

  // f_d_:
  f_h_t f_h_;
  f_d_t f_d_;

  // signals act always on velocity grid points.
  // since we can have different acting locations, we need to have a 1dview
  // to store the VpGlobal ID of all grid points where sources act on
  vp_gids_h_t myVpGids_h_;
  vp_gids_d_t myVpGids_d_;

  // max freq of all sources
  sc_t maxFreq_ = {};
  // time step size
  const sc_t dt_ = {};

public:
  template <class parser_t, class mesh_info_t, class app_t, typename param_t>
  RankTwoForcing(const signal_instances_h_type & signals,
		 const parser_t & parser,
		 const mesh_info_t & meshInfo,
		 const app_t & appObj,
		 const param_t & depthsKm,
		 const param_t & anglesDeg,
		 bool scaleByDt = false)
    : signals_(signals),
      f_h_("Fh", signals.extent(0)),
      f_d_("Fd", signals.extent(0)),
      myVpGids_h_("vpGidsH", signals.extent(0)),
      myVpGids_d_("vpGidsD", signals.extent(0)),
      dt_(parser.getTimeStepSize())
  {
    KokkosBlas::fill(f_h_, constants<sc_t>::zero());
    KokkosBlas::fill(f_d_, constants<sc_t>::zero());

    const auto gidsVp = appObj.viewGidListHost(dofId::vp);
    const auto coords = appObj.viewCoordsHost(dofId::vp);

    // where this forcing is located
    const auto domainSurfaceRadiusKm = meshInfo.getMaxRadiusKm();
    for (std::size_t i=0; i<depthsKm.size(); ++i)
    {
      const auto myRadiusKm = domainSurfaceRadiusKm - depthsKm[i];
      const auto myAngleDeg = anglesDeg[i];

      // find the vpGid identifying the grid point where the source is mapped to
      mapPointSourceToGridPoint(myAngleDeg, myRadiusKm, depthsKm[i],
				meshInfo.viewDomainBounds(),
				meshInfo.getNumVpPts(), gidsVp, coords,
				meshInfo.getAngularSpacing(),
				myVpGids_h_(i));
    }
    Kokkos::deep_copy(myVpGids_d_, myVpGids_h_);
    computeMaxFrequency(signals);

    // for (std::size_t i=0; i<signals.extent(0); ++i)
    // {
    //   const auto & signalIt = signals(i);
    //   auto time = constants<sc_t>::zero();
    //   for (std::size_t iStep = 1; iStep<=f_h_.extent(0); ++iStep){
    // 	signalIt(time, f_h_(iStep-1, i));
    // 	time = iStep * dt_;
    //   }
    // }
  }

  // for single source, max frquency is the frequency of the source signal
  sc_t getMaxFreq() const{ return maxFreq_; }

  vp_gids_h_t getVpGidsHost() const{ return myVpGids_h_; }
  vp_gids_d_t getVpGidsDevice() const{ return myVpGids_d_; }

  f_d_t viewForcingDevice() const{ return f_d_; }

  auto getForcingAtStep(const std::size_t & step) const{
    return Kokkos::subview(f_h_, step-1, Kokkos::ALL());
  }

  void evaluate(const sc_t & time, const std::size_t & step)
  {
    for (std::size_t i=0; i<signals_.extent(0); ++i)
    {
      const auto & signalIt = signals_(i);
      signalIt(time, f_h_(i));
    }
    Kokkos::deep_copy(f_d_, f_h_);

    // // KokkosBlas::fill(f_d_, constants<sc_t>::zero());
    // // for (std::size_t j=0; j<f_h_.extent(1); ++j){
    // const auto src = Kokkos::subview(f_h_, step-1, Kokkos::ALL());
    // //const auto des = Kokkos::subview(f_d_, myVpGids_(j), j);
    // Kokkos::deep_copy(f_d_, src);
    // // }
  }

  void complexityOfEvaluateMethod(double & memCostMB, double & flopsCost) const
  {
    // no operation is done during evaluate, just copying, see above
    const double memMBCostFill=0.;
    //1.*( f_d_.extent(0)*f_d_.extent(1)*sizeof(sc_t) )/1024./1024.;
    const double memMBCostCopy=0.;
    //.1.*( f_d_.extent(1)*sizeof(sc_t) )/1024./1024.;
    memCostMB = memMBCostFill + memMBCostCopy;
    flopsCost = 0.;
  }

private:
  template<class signals_t>
  void computeMaxFrequency(const signals_t & signals)
  {
    maxFreq_ = std::numeric_limits<sc_t>::min();
    for (std::size_t i=0; i<signals.extent(0); ++i){
      maxFreq_ = std::max( maxFreq_, signals[i].getFrequency() );
    }
  }
};

#endif
