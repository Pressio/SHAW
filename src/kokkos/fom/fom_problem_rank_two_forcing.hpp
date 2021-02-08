
#ifndef FOM_PROBLEM_RANK_TWO_FORCING_HPP_
#define FOM_PROBLEM_RANK_TWO_FORCING_HPP_

namespace kokkosapp{

class FomProblemRankTwoForcing
  : public kokkosapp::commonTypes
{
  using kokkosapp::commonTypes::scalar_t;
  using kokkosapp::commonTypes::sc_t;
  using kokkosapp::commonTypes::int_t;
  using kokkosapp::commonTypes::parser_t;
  using kokkosapp::commonTypes::mesh_info_t;
  using kokkosapp::commonTypes::klr;
  using kokkosapp::commonTypes::kll;
  using kokkosapp::commonTypes::exe_space;

  static constexpr bool usingFullMesh = true;
  using state_d_t	= Kokkos::View<sc_t**, klr, exe_space>;
  using state_h_t	= typename state_d_t::host_mirror_type;
  using jacobian_d_type = KokkosSparse::CrsMatrix<sc_t, int_t, exe_space>;
  using fom_t		= ShWavePP<sc_t, int_t, mesh_info_t, jacobian_d_type, exe_space>;
  using obs_t		= StateObserver<int_t, sc_t>;
  using seismogram_t	= Seismogram<int_t, sc_t>;
  using forcing_t       = RankTwoForcing<sc_t, state_d_t, int_t>;

private:
  // parser with inputs
  const parser_t & parser_;

  // object with info about the mesh
  const mesh_info_t & meshInfo_;

  // material model object
  std::shared_ptr<MaterialModelBase<scalar_t>> materialObj_;

  const int_t fSize_;
  const int_t nVp_;
  const int_t nSp_;

  fom_t appObj_;
  state_d_t xVp_d_;
  state_d_t xSp_d_;
  obs_t observerObj_;

public:
  FomProblemRankTwoForcing(const parser_t & parser,
			   const mesh_info_t & meshInfo,
			   std::shared_ptr<MaterialModelBase<scalar_t>> materialObj)
    : parser_(parser),
      meshInfo_(meshInfo),
      materialObj_(materialObj),
      fSize_(parser.getForcingSize()),
      nVp_(meshInfo_.getNumVpPts()),
      nSp_(meshInfo_.getNumSpPts()),
      appObj_(meshInfo_),
      xVp_d_("xVp_d", nVp_, fSize_),
      xSp_d_("xSp_d", nSp_, fSize_),
      observerObj_(nVp_, nSp_, parser, fSize_)
  {}

public:
  void execute()
  {
    multiRunSamplingForcingPeriod();
  }

private:
  void multiRunSamplingForcingPeriod()
  {
    /* here we sample the forcing period,
       which means that:
       1. the material does not change
       2. the other properties (like location) of the source do not change
    */

    appObj_.computeJacobians(*materialObj_);

    // seismogram
    seismogram_t seismoObj(parser_, meshInfo_, appObj_, fSize_);

    // create vector of signals
    const auto periods = parser_.getValues(0);
    std::vector<Signal<sc_t>> signals;
    for (auto i=0; i<periods.size(); ++i){
      // use the signal that was set from input file
      // so everything remains the same except for the
      signals.emplace_back(parser_.getSignal());
      signals.back().resetPeriod(periods[i]);
    }

    std::cout << "Doing FOM with sampling of forcing period" << std::endl;
    std::cout << "Total number of samples " << signals.size() << std::endl;

    // create a forcing object, this does mem allocation
    // (in loop below, only thing that changes is the signal NOT the location,
    // so it is fine to create the nominal forcing and the in the loop replace signal)
    forcing_t forcing(parser_, meshInfo_, appObj_);

    const std::size_t numSets = signals.size()/fSize_;
    for (std::size_t i=0; i<numSets; ++i)
    {
      // each loop iteration handles fSize_ signals
      // so replace signals starting from i*fSize_
      forcing.replaceSignals(signals, i*fSize_);

      // need to recheck that the new signal still meets conditions
      doChecks(forcing);

      // reset observer and seismogram
      observerObj_.prepForNewRun(i);
      seismoObj.prepForNewRun(i);

      // run fom
      runFom(parser_.getNumSteps(), parser_.getTimeStepSize(),
      	     appObj_, forcing, observerObj_, seismoObj, xVp_d_, xSp_d_);

      processCollectedData(seismoObj);
    }

    // coordinates only need to be written once
    processCoordinates();
  }

  template <typename forcing_t>
  void doChecks(const forcing_t & forcing){
    if (parser_.checkDispersion())
      checkDispersionCriterion(meshInfo_, forcing.getMaxFreq(),
    			       appObj_.getMinShearWaveVelocity());

    if (parser_.checkCfl())
      checkCfl(meshInfo_, parser_.getTimeStepSize(), appObj_.getMaxShearWaveVelocity());
  }

  void processCoordinates()
  {
    if(parser_.enableSnapshotMatrix()){
      appObj_.writeCoordinatesToFile(dofId::vp);
      appObj_.writeCoordinatesToFile(dofId::sp);
    }
  }

  template <typename seismo_t>
  void processCollectedData(const seismo_t & seismoObj)
  {
    const auto startTime  = std::chrono::high_resolution_clock::now();

    if(parser_.enableSnapshotMatrix()){
      observerObj_.writeSnapshotMatrixToFile(dofId::vp);
      observerObj_.writeSnapshotMatrixToFile(dofId::sp);
    }

    if(parser_.enableSeismogram()){
      seismoObj.writeReceiversToFile();
    }

    const auto finishTime = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double> elapsed = finishTime - startTime;
    std::cout << "\nfinalProcessTime = " << std::fixed << std::setprecision(10)
	      << elapsed.count();
    std::cout << "\n";
  }

};

}//end namespace
#endif
