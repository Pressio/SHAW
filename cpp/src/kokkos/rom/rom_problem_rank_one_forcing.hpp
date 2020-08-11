
#ifndef ROM_PROBLEM_RANK_ONE_FORCING_HPP_
#define ROM_PROBLEM_RANK_ONE_FORCING_HPP_

#include "../../shared/all.hpp"
#include "../shwavepp.hpp"
#include "../common_types.hpp"
#include "compute_rom_jacobians.hpp"
#include "load_basis.hpp"
#include "run_rom_rank_one_forcing.hpp"
#include "rom_problem_base.hpp"

namespace kokkosapp{

class RomProblemRankOneForcing final
  : public RomProblemBase, kokkosapp::commonTypes
{
  using kokkosapp::commonTypes::scalar_t;
  using kokkosapp::commonTypes::sc_t;
  using kokkosapp::commonTypes::int_t;
  using kokkosapp::commonTypes::parser_t;
  using kokkosapp::commonTypes::klr;
  using kokkosapp::commonTypes::kll;
  using kokkosapp::commonTypes::exe_space;

  static constexpr bool usingFullMesh = true;
  using mesh_info_t	= MeshInfo<sc_t, int_t, usingFullMesh>;

  // here the state is a rank-1 view so layout should not matter
  using state_d_t	= Kokkos::View<sc_t*, kll, exe_space>;
  using state_h_t	= typename state_d_t::host_mirror_type;
  using jacobian_d_type = KokkosSparse::CrsMatrix<sc_t, int_t, exe_space>;
  using fom_t		= ShWavePP<sc_t, int_t, mesh_info_t, jacobian_d_type, exe_space>;

  // here the basis is a rank-2 view and layout matters for binary IO
  // but not for performance, the basis are only need to compute reduced jacobians
  // leave it to kll which is kokkos::layoutLeft
  using basis_d_t      = Kokkos::View<sc_t**, kll, exe_space>;
  using basis_h_t      = typename basis_d_t::host_mirror_type;

  // here the rom_jac is a rank-2 view and layout matters for performance
  // should pick ll/lr to optimize gemv for every time step
  using rom_jac_d_t    = Kokkos::View<sc_t**, klr, exe_space>;
  using rom_jac_h_t    = typename rom_jac_d_t::host_mirror_type;

  using obs_t		= Observer<int_t, sc_t, state_d_t>;
  using forcing_t       = RankOneForcing<sc_t, state_d_t, int_t>;

private:
  const parser_t & parser_;
  const mesh_info_t meshInfo_;
  const int_t nVpFom_;
  const int_t nSpFom_;
  const int_t nVp_;
  const int_t nSp_;

  state_d_t xVp_d_;
  state_d_t xSp_d_;
  basis_d_t phiVp_d_;
  basis_d_t phiSp_d_;
  rom_jac_d_t Jvp_d_;
  rom_jac_d_t Jsp_d_;
  obs_t observerObj_;

public:
  RomProblemRankOneForcing(const parser_t & parser)
    : parser_(parser),
      meshInfo_(parser.getMeshDir()),
      nVpFom_(meshInfo_.getNumVpPts()),
      nSpFom_(meshInfo_.getNumSpPts()),
      nVp_(parser.getRomSize(dofId::vp)),
      nSp_(parser.getRomSize(dofId::sp)),
      xVp_d_("xVp_d", nVp_),
      xSp_d_("xSp_d", nSp_),
      phiVp_d_("phiVp_d", nVpFom_, nVp_),
      phiSp_d_("phiSp_d", nSpFom_, nSp_),
      Jvp_d_("romJvp_d", nVp_, nSp_),
      Jsp_d_("romJsp_d", nSp_, nVp_),
      observerObj_(nVp_, nSp_, parser)
  {
    loadBasis<int_t, sc_t, basis_h_t>(parser_, phiVp_d_, phiSp_d_);
  }

public:
  void execute() final
  {
    if (parser_.enableSampling())
    {
      const auto param = parser_.getNameParamToSample(0);
      if (param == samplable::signalPeriod)
	multiRunSamplingForcingPeriod();
      else{
	const auto msg = "rom:rank1forcing: sampling for param!=signalPeriod not yet supported";
	throw std::runtime_error(msg);
      }
    }
    else{
      executeSingleRun();
    }
  }

private:
  void executeSingleRun()
  {
    // compute rom Jacobians
    auto matObj = createMaterialModel<sc_t>(parser_);
    fom_t fomObj(meshInfo_);
    fomObj.computeJacobiansWithMatProp(*matObj);
    computeRomOperatorsUsingFomJacsWithMatProp(parser_, fomObj, phiVp_d_, phiSp_d_, Jvp_d_, Jsp_d_);

    // forcing object
    forcing_t forcingObj(parser_, meshInfo_, fomObj);

    // index of where the source acts
    const auto fIndex   = forcingObj.getVpGid();

    // create vector containg phiVp^T*rhoInv which pre-multiplies the forcing
    state_d_t phiVpRhoInvVec_d("phiVpRhoInv_d", nVp_);
    // extract the target row from phiVp
    auto phiVpRow_d = Kokkos::subview(phiVp_d_, fIndex, Kokkos::ALL());

    // copy phiVp(fIndex,:) to phiVpRhoInvVec
    Kokkos::deep_copy(phiVpRhoInvVec_d, phiVpRow_d);

    // extract the target value from rhoInv
    const auto rhoInvVpValue = fomObj.viewInvDensityHost(dofId::vp)(fIndex);
    // scale
    KokkosBlas::scal(phiVpRhoInvVec_d, rhoInvVpValue, phiVpRhoInvVec_d);

    // need to recheck that the new signal still meets conditions
    doChecks(forcingObj, fomObj);

    kokkosapp::runRomRankOneForcing(parser_.getNumSteps(), parser_.getTimeStepSize(),
    				    phiVpRhoInvVec_d, forcingObj, observerObj_,
    				    Jvp_d_, Jsp_d_, xVp_d_, xSp_d_);
    processCollectedData();
  }

  void multiRunSamplingForcingPeriod()
  {
    /* here we sample the forcing period,
       which means that:
       1. the material does not change
       2. the other properties (like location) of the source do not change
    */

    /*
     * create and store material prop
     * only do it once since material does not change
    */

    auto matObj = createMaterialModel<sc_t>(parser_);
    // compute rom Jacobians: only onces, these don't change
    fom_t fomObj(meshInfo_);
    fomObj.computeJacobiansWithMatProp(*matObj);
    computeRomOperatorsUsingFomJacsWithMatProp(parser_, fomObj, phiVp_d_, phiSp_d_, Jvp_d_, Jsp_d_);

    // create vector of signals
    const auto periods = parser_.getValues(0);
    using signal_t = Signal<sc_t>;
    std::vector<signal_t> signals;
    for (auto i=0; i<periods.size(); ++i)
    {
      signals.emplace_back(parser_.getSignal());
      signals.back().resetPeriod(periods[i]);
    }

    std::cout << "Doing ROM with sampling of forcing period" << std::endl;
    std::cout << "Total number of samples " << signals.size() << std::endl;

    // create a forcing object with mem allocation
    // (in loop below, only thing that changes is
    // the signal NOT the location of the signal, so it is fine to
    // create the nominal forcing and the in the loop below replace signal)
    forcing_t forcingObj(parser_, meshInfo_, fomObj);

    // index of where the source acts
    const auto fIndex   = forcingObj.getVpGid();

    // create vector containg phiVp^T*rhoInv which pre-multiplies the forcing
    state_d_t phiVpRhoInvVec_d("phiVpRhoInv_d", nVp_);
    // extract the target row from phiVp
    auto phiVpRow_d = Kokkos::subview(phiVp_d_, fIndex, Kokkos::ALL());

    // copy phiVp(fIndex,:) to phiVpRhoInvVec
    Kokkos::deep_copy(phiVpRhoInvVec_d, phiVpRow_d);

    // extract the target value from rhoInv
    const auto rhoInvVpValue = fomObj.viewInvDensityHost(dofId::vp)(fIndex);
    // scale
    KokkosBlas::scal(phiVpRhoInvVec_d, rhoInvVpValue, phiVpRhoInvVec_d);

    // loop over signals
    for (std::size_t iSample=0; iSample<signals.size(); ++iSample)
    {
      // replace signal (no new allocations happen here)
      forcingObj.replaceSignal(signals[iSample]);

      // need to recheck that the new signal still meets conditions
      doChecks(forcingObj, fomObj);

      // reset observer
      observerObj_.prepForNewRun(iSample);

      // run
      kokkosapp::runRomRankOneForcing(parser_.getNumSteps(), parser_.getTimeStepSize(),
				      phiVpRhoInvVec_d, forcingObj, observerObj_,
				      Jvp_d_, Jsp_d_, xVp_d_, xSp_d_);
      processCollectedData(iSample);
    }
  }

  template <typename forcing_t, typename fom_t>
  void doChecks(const forcing_t & forcing, fom_t & fomObj){
    if (parser_.checkDispersion())
      checkDispersionCriterion(meshInfo_, forcing.getMaxFreq(),
    			       fomObj.getMinShearWaveVelocity());

    if (parser_.checkCfl())
      checkCfl(meshInfo_, parser_.getTimeStepSize(), fomObj.getMaxShearWaveVelocity());
  }

  void processCollectedData(std::size_t iSample = 0)
  {
    // only process final/stored data if dummy basis = false
    if (parser_.enableRandomDummyBasis() == false)
    {

      const auto startTime  = std::chrono::high_resolution_clock::now();

      if(parser_.enableSnapshotMatrix()){
	observerObj_.writeSnapshotMatrixToFile(dofId::vp);
	observerObj_.writeSnapshotMatrixToFile(dofId::sp);
      }

      const auto finishTime = std::chrono::high_resolution_clock::now();
      const std::chrono::duration<double> elapsed = finishTime - startTime;
      std::cout << "\nfinalProcessTime = " << std::fixed << std::setprecision(10)
		<< elapsed.count();
      std::cout << "\n";

      // if(parser_.enableWriteFinalState()){
      // 	processFinalStateToFile(dofId::vp, iSample);
      // 	processFinalStateToFile(dofId::sp, iSample);
      // }
    }
  }

  // void processFinalStateToFile(const dofId dof, std::size_t iSample = 0)
  // {
  //   constexpr auto one	= constants<sc_t>::one();
  //   constexpr auto zero	= constants<sc_t>::zero();

  //   const auto dofName	  = dofIdToString(dof);
  //   const auto romState_d = (dof==dofId::vp) ? xVp_d_ : xSp_d_;
  //   const auto phi_d	  = (dof==dofId::vp) ? phiVp_d_ : phiSp_d_;

  //   const auto fomStateSize     = (dof==dofId::vp) ? nVpFom_ : nSpFom_;
  //   const auto fomStateFileName = "finalFomState_"+dofName+"_"+std::to_string(iSample);
  //   const auto romStateFileName = "finalRomState_"+dofName+"_"+std::to_string(iSample);

  //   state_d_t fomState_d("fomState_d", fomStateSize);
  //   const char ct_N	= 'N';
  //   KokkosBlas::gemv(&ct_N, one, phi_d, romState_d, zero, fomState_d);

  //   auto fomState_h = Kokkos::create_mirror_view(fomState_d);
  //   auto romState_h = Kokkos::create_mirror_view(romState_d);
  //   Kokkos::deep_copy(fomState_h, fomState_d);
  //   Kokkos::deep_copy(romState_h, romState_d);

  //   writeToFile(fomStateFileName, fomState_h, parser_.writeFinalStateBinary());
  //   // the rom state is small, so use ascii
  //   writeToFile(romStateFileName, romState_h, false);
  // }

};

}//end namespace
#endif
