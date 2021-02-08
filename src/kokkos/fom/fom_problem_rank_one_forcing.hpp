
#ifndef FOM_PROBLEM_RANK_ONE_FORCING_HPP_
#define FOM_PROBLEM_RANK_ONE_FORCING_HPP_

namespace kokkosapp{

template<typename T>
class FomProblemRankOneForcing
{
public:
  using scalar_type     = typename T::scalar_type;
  using parser_type     = typename T::parser_type;
  using mesh_info_type  = typename T::mesh_info_type;
  using state_d_type    = typename T::state_d_type;
  using forcing_type    = typename T::forcing_type;
  using observer_type   = typename T::observer_type;
  using seismogram_type = typename T::seismogram_type;
  using mesh_ord_type	= typename mesh_info_type::ordinal_type;

private:
  // parser with inputs
  const parser_type & parser_;

  // object with info about the mesh
  const mesh_info_type & meshInfo_;

  // material model object
  std::shared_ptr<MaterialModelBase<scalar_type>> materialObj_;

  // number of velocity DOFs
  const mesh_ord_type nVp_;
  // number of stresses DOFs
  const mesh_ord_type nSp_;

  // system object to construct operators
  ShWavePP<T> appObj_;
  // state vector for the velocity DOFs
  state_d_type xVp_d_;
  // state vector for the stress DOFs
  state_d_type xSp_d_;
  // observer object to monitor the time evolution
  observer_type observerObj_;

public:
  FomProblemRankOneForcing(const parser_type & parser,
			   const mesh_info_type & meshInfo,
			   std::shared_ptr<MaterialModelBase<scalar_type>> materialObj)
    : parser_(parser),
      meshInfo_(meshInfo),
      materialObj_(materialObj),
      nVp_(meshInfo_.getNumVpPts()),
      nSp_(meshInfo_.getNumSpPts()),
      appObj_(meshInfo_),
      xVp_d_("xVp_d", nVp_),
      xSp_d_("xSp_d", nSp_),
      observerObj_(nVp_, nSp_, parser)
  {}

public:
  void execute()
  {
    // if (parser_.enableSampling())
    // {
    //   const auto param = parser_.getNameParamToSample(0);
    //   if (param == samplable::signalPeriod)
    // 	multiRunSamplingForcingPeriod();
    //   else{
    // 	const auto msg = "fom:rank1forcing: sampling for param!=signalPeriod not yet supported";
    // 	throw std::runtime_error(msg);
    //   }
    // }
    // else{
    singleRun();
      //}
  }

private:
  void singleRun()
  {
    // use material model to compute Jacobian matrices
    appObj_.computeJacobians(*materialObj_);

    // seismogram: stores the seismogram at locations specified in input file
    seismogram_type seismoObj(parser_, meshInfo_, appObj_);

    // construct forcing using signal info from parser
    forcing_type forcing(parser_, meshInfo_, appObj_);

    // run checks
    doChecks(forcing);

    // run fom
    runFom(parser_.getNumSteps(), parser_.getTimeStepSize(),
	   appObj_, forcing, observerObj_, seismoObj,
	   xVp_d_, xSp_d_);

    processCoordinates();
    processCollectedData(seismoObj);
  }

  // void multiRunSamplingForcingPeriod()
  // {
  //   /* here we sample the forcing period,
  //      which means that:
  //      1. the material does not change
  //      2. the other properties (like location) of the source do not change
  //   */

  //   /*
  //    * create and store material prop
  //    * only do it once since material does not change
  //   */
  //   auto matObj = createMaterialModel<scalar_type>(parser_, meshInfo_);
  //   appObj_.computeJacobiansWithMatProp(*matObj);

  //   // seismogram
  //   seismogram_type seismoObj(parser_, meshInfo_, appObj_);

  //   // create vector of signals
  //   const auto periods = parser_.getValues(0);
  //   using signal_t = Signal<scalar_type>;
  //   std::vector<signal_t> signals;
  //   for (auto i=0; i<periods.size(); ++i)
  //   {
  //     signals.emplace_back(parser_.getSignal());
  //     signals.back().resetPeriod(periods[i]);
  //   }

  //   std::cout << "Doing FOM with sampling of forcing period" << std::endl;
  //   std::cout << "Total number of samples " << signals.size() << std::endl;

  //   // create a forcing object with mem allocation
  //   // (in loop below, only thing that changes is
  //   // the signal NOT the location of the signal, so it is fine to
  //   // create the nominal forcing and the in the loop below replace signal)
  //   forcing_type forcing(parser_, meshInfo_, appObj_);

  //   // loop over signals
  //   for (std::size_t iSample=0; iSample<signals.size(); ++iSample)
  //   {
  //     // replace signal (no new allocations happen here)
  //     forcing.replaceSignal(signals[iSample]);

  //     // need to recheck that the new signal still meets conditions
  //     doChecks(forcing);

  //     // reset observer and seismogram
  //     observerObj_.prepForNewRun(iSample);

  //     // run fom
  //     runFom(true, parser_.exploitForcingSparsity(),
  // 	     parser_.getNumSteps(), parser_.getTimeStepSize(),
  // 	     appObj_, forcing, observerObj_, seismoObj, xVp_d_, xSp_d_);

  //     processCollectedData(seismoObj, iSample);
  //   }

  //   // coordinates only need to be written once
  //   processCoordinates();
  // }

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
  void processCollectedData(const seismo_t & seismoObj, std::size_t iSample = 0)
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
    std::cout << "\nfinalProcessTime = " << std::fixed << std::setprecision(10) << elapsed.count();
    std::cout << "\n";
  }
};

}//end namespace
#endif
