
#ifndef FOM_PROBLEM_RANK_ONE_FORCING_HPP_
#define FOM_PROBLEM_RANK_ONE_FORCING_HPP_

#include "../../shared/all.hpp"
#include "../shwavepp.hpp"
#include "../common_types.hpp"
#include "run_fom.hpp"

namespace kokkosapp{

class FomProblemRankOneForcing
  : public kokkosapp::commonTypes
{
public:
  using kokkosapp::commonTypes::scalar_t;
  using kokkosapp::commonTypes::sc_t;
  using kokkosapp::commonTypes::int_t;
  using kokkosapp::commonTypes::parser_t;
  using kokkosapp::commonTypes::klr;
  using kokkosapp::commonTypes::kll;
  using kokkosapp::commonTypes::exe_space;

  using state_d_t	= Kokkos::View<sc_t*, kll, exe_space>;
  using state_h_t	= typename state_d_t::host_mirror_type;
  using jacobian_d_type = KokkosSparse::CrsMatrix<sc_t, int_t, exe_space>;
  using mesh_info_t	= MeshInfo<sc_t, int_t>;
  using fom_t		= ShWavePP<sc_t, int_t, mesh_info_t, jacobian_d_type, exe_space>;
  using forcing_t       = RankOneForcing<sc_t, state_d_t, int_t>;
  using obs_t		= StateObserver<int_t, sc_t>;
  using seismogram_t	= Seismogram<int_t, sc_t>;

private:
  // parser with inputs
  const parser_t & parser_;

  // object with info about the mesh
  const mesh_info_t meshInfo_;

  // number of velocity DOFs
  const int_t nVp_;
  // number of stresses DOFs
  const int_t nSp_;

  // system object to construct operators
  fom_t appObj_;
  // state vector for the velocity DOFs
  state_d_t xVp_d_;
  // state vector for the stress DOFs
  state_d_t xSp_d_;
  // observer object to monitor the time evolution
  obs_t observerObj_;

public:
  FomProblemRankOneForcing(const parser_t & parser)
    : parser_(parser),
      meshInfo_(parser.getMeshDir()),
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
    // create material model object
    auto matObj = createMaterialModel<sc_t>(parser_);

    // use material model to compute Jacobian matrices
    appObj_.computeJacobians(*matObj);

    // seismogram: stores the seismogram at locations specified in input file
    seismogram_t seismoObj(parser_, meshInfo_, appObj_);

    // construct forcing using signal info from parser
    forcing_t forcing(parser_, meshInfo_, appObj_);

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
  //   auto matObj = createMaterialModel<sc_t>(parser_);
  //   appObj_.computeJacobiansWithMatProp(*matObj);

  //   // seismogram
  //   seismogram_t seismoObj(parser_, meshInfo_, appObj_);

  //   // create vector of signals
  //   const auto periods = parser_.getValues(0);
  //   using signal_t = Signal<sc_t>;
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
  //   forcing_t forcing(parser_, meshInfo_, appObj_);

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
