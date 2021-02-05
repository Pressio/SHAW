
#ifndef RANK_TWO_FORCING_HPP_
#define RANK_TWO_FORCING_HPP_

#include "map_point_source_to_velocity_grid_point.hpp"
#include "KokkosBlas1_fill.hpp"
#include "KokkosBlas1_scal.hpp"

template <typename sc_t, typename state_d_t, typename int_t>
class RankTwoForcing
{
  static_assert(is_kokkos_2dview< state_d_t >::value,
		"Rank-2 forcing must use a rank-2 kokkos view");

  using state_h_t = typename state_d_t::host_mirror_type;

  // f_h_ is a 2d view where
  // each column contains the full time-series for one signal instance
  // therefore: rows index the time, cols refer to different signal instances
  state_h_t f_h_;

  // f_d_: is a 2d view, where each col represents the full full forcing vector
  // at a given time, over the full mesh (so f_d_ can be large)
  // Note that: f_d_ must be updated at every time step
  state_d_t f_d_;

  // A signal acts always on a velocity point.
  // myVpGid_ identifies which velocity point all signals act on.
  int_t myVpGid_ = -1;

  sc_t maxFreq_ = {};
  const sc_t dt_ = {};
  const int_t NSteps_ = {};
  bool scaleByDt_ = false;

public:
  template <typename parser_t, typename mesh_info_t, typename app_t>
  RankTwoForcing(const parser_t & parser,
		 const mesh_info_t & meshInfo,
		 const app_t & appObj,
		 const sc_t depthKm,
		 const sc_t angleDeg,
		 bool scaleByDt = false)
    : f_h_("Fh", parser.getNumSteps(), parser.getForcingSize()),
      f_d_("Fd", meshInfo.getNumVpPts(), parser.getForcingSize()),
      dt_(parser.getTimeStepSize()),
      NSteps_(parser.getNumSteps()),
      scaleByDt_(scaleByDt)
  {
    const auto gidsVp = appObj.viewGidListHost(dofId::vp);
    const auto coords = appObj.viewCoordsHost(dofId::vp);

    // where this forcing is located
    const auto domainSurfaceRadiusKm = meshInfo.getMaxRadiusKm();
    const sc_t myRadiusKm = domainSurfaceRadiusKm - depthKm;

    // find the vpGid identifying the grid point where the source is mapped to
    mapPointSourceToGridPoint(angleDeg, myRadiusKm, depthKm,
  			      meshInfo.viewDomainBounds(),
  			      meshInfo.getNumVpPts(), gidsVp, coords,
  			      meshInfo.getAngularSpacing(), myVpGid_);
  }

  template <typename parser_t, typename mesh_info_t, typename app_t>
  RankTwoForcing(const parser_t    & parser,
		 const mesh_info_t & meshInfo,
		 const app_t	   & appObj)
    : RankTwoForcing(parser, meshInfo, appObj,
		     parser.getSourceProperty("depth"),
		     parser.getSourceProperty("angle"))
  {}

  template <typename signal_t>
  void replaceSignals(const std::vector<signal_t> & signals,
		      std::size_t startIndex)
  {
    // maxFrequency
    maxFreq_ = std::numeric_limits<sc_t>::min();

    assert(signals.size() == f_h_.extent(1));
    KokkosBlas::fill(f_h_, constants<sc_t>::zero());

    std::size_t c=0;
    for (std::size_t i=startIndex; i<startIndex+f_h_.extent(1); ++i)
    {
      const auto & signalIt = signals[i];

      maxFreq_ = std::max( maxFreq_, signalIt.getFrequency() );

      // store the full time series of the signal into the host array
      sc_t time = constants<sc_t>::zero();
      for (int_t iStep = 1; iStep<=NSteps_; ++iStep){
  	signalIt(time, f_h_(iStep-1, c));
  	time = iStep * dt_;
      }
      c++;
    }
    if (scaleByDt_)
      KokkosBlas::scal(f_h_, dt_, f_h_);
  }

public:
  // for single source, max frquency is the frequency of the source signal
  sc_t getMaxFreq() const{ return maxFreq_; }

  int_t getVpGid() const{ return myVpGid_; }

  state_d_t viewForcingDevice() const{ return f_d_; }

  auto getForcingAtStep(std::size_t step) const{
    return Kokkos::subview(f_h_, step-1, Kokkos::ALL());
  }

  void evaluate(sc_t time, std::size_t step){
    KokkosBlas::fill(f_d_, constants<sc_t>::zero());
    const auto src = Kokkos::subview(f_h_, step-1, Kokkos::ALL());
    const auto des = Kokkos::subview(f_d_, myVpGid_, Kokkos::ALL());
    Kokkos::deep_copy(des, src);
  }

  void complexityOfEvaluateMethod(double & memCostMB, double flopsCost) const{
    // no operation is done during evaluate, just copying, see above
    const double memMBCostFill = 1.*( f_d_.extent(0)*f_d_.extent(1)*sizeof(sc_t) )/1024./1024.;
    const double memMBCostCopy = 1.*( f_d_.extent(1)*sizeof(sc_t) )/1024./1024.;
    memCostMB = memMBCostFill + memMBCostCopy;
    flopsCost = 0.;
  }
};

#endif
