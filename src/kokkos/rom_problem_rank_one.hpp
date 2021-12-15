/*
//@HEADER
// ************************************************************************
//
// rom_problem_rank_one_forcing.hpp
//                     		Pressio/SHAW
//                         Copyright 2019
// National Technology & Engineering Solutions of Sandia, LLC (NTESS)
//
// Under the terms of Contract DE-NA0003525 with NTESS, the
// U.S. Government retains certain rights in this software.
//
// Pressio is licensed under BSD-3-Clause terms of use:
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
// contributors may be used to endorse or promote products derived
// from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Francesco Rizzi (fnrizzi@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef ROM_PROBLEM_RANK_ONE_FORCING_HPP_
#define ROM_PROBLEM_RANK_ONE_FORCING_HPP_

namespace kokkosapp{

template<typename T>
class RomProblemRankOneForcing
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

  using basis_d_type   = typename T::basis_d_type;
  using basis_h_type   = typename T::basis_h_type;
  using rom_jac_d_type = typename T::rom_jac_d_type;
  using rom_jac_h_type = typename T::rom_jac_h_type;

private:
  const parser_type & parser_;
  const mesh_info_type meshInfo_;
  const mesh_ord_type nVpFom_;
  const mesh_ord_type nSpFom_;
  const mesh_ord_type nVp_;
  const mesh_ord_type nSp_;
  ShWavePP<T> appObj_;
  state_d_type xVp_d_;
  state_d_type xSp_d_;
  observer_type observerObj_;

  basis_d_type phiVp_d_;
  basis_d_type phiSp_d_;
  rom_jac_d_type Jvp_d_;
  rom_jac_d_type Jsp_d_;

public:
  RomProblemRankOneForcing(const parser_type & parser,
			   const mesh_info_type & meshInfo,
			   const MaterialModelBase<scalar_type> & materialObj)
    : parser_(parser),
      meshInfo_(meshInfo),
      nVpFom_(meshInfo_.getNumVpPts()),
      nSpFom_(meshInfo_.getNumSpPts()),
      nVp_(parser.getRomSize(dofId::vp)),
      nSp_(parser.getRomSize(dofId::sp)),
      /* create app object */
      appObj_(meshInfo_, materialObj),
      /* states */
      xVp_d_("xVp_d", nVp_),
      xSp_d_("xSp_d", nSp_),
      /* observer */
      observerObj_(nVp_, nSp_, parser),
      /* basis */
      phiVp_d_("phiVp_d", nVpFom_, nVp_),
      phiSp_d_("phiSp_d", nSpFom_, nSp_),
      /* rom jacobians */
      Jvp_d_("romJvp_d", nVp_, nSp_),
      Jsp_d_("romJsp_d", nSp_, nVp_)
  {
    fillBasis<mesh_ord_type, scalar_type>(parser_, phiVp_d_, phiSp_d_);
  }

public:
  void operator()()
  {
    executeSingleRun();
  }

private:
  void executeSingleRun()
  {
    // compute rom Jacobians
    computeRomOperatorsUsingFomJacsWithMatProp(parser_, appObj_,
					       phiVp_d_, phiSp_d_,
					       Jvp_d_, Jsp_d_);

    // forcing object
    forcing_type forcingObj(parser_, meshInfo_, appObj_);

    // index of where the source acts
    const auto fIndex = forcingObj.getVpGid();

    // create vector containg phiVp^T*rhoInv which pre-multiplies the forcing
    state_d_type phiVpRhoInvVec_d("phiVpRhoInv_d", nVp_);
    // extract the target row from phiVp
    auto phiVpRow_d = Kokkos::subview(phiVp_d_, fIndex, Kokkos::ALL());

    // copy phiVp(fIndex,:) to phiVpRhoInvVec
    Kokkos::deep_copy(phiVpRhoInvVec_d, phiVpRow_d);

    // extract the target value from rhoInv
    const auto rhoInvVpValue = appObj_.viewInvDensityHost(dofId::vp)(fIndex);
    // scale
    KokkosBlas::scal(phiVpRhoInvVec_d, rhoInvVpValue, phiVpRhoInvVec_d);

    // run checks
    checkCflCondition();
    checkDispersion(forcingObj.getMaxFreq());

    runRomRankOneForcing(parser_.getNumSteps(),
			 parser_.getTimeStepSize(),
			 phiVpRhoInvVec_d, forcingObj,
			 observerObj_,
			 Jvp_d_, Jsp_d_, xVp_d_, xSp_d_);
    processCollectedData();
  }

  void checkDispersion(const scalar_type & freq){
    if (parser_.checkDispersion()){
      checkDispersionCriterion(meshInfo_, freq, appObj_.getMinShearWaveVelocity());
    }
  }

  void checkCflCondition(){
    if (parser_.checkCfl()){
      checkCfl(meshInfo_, parser_.getTimeStepSize(), appObj_.getMaxShearWaveVelocity());
    }
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
      std::cout << "\nfinalProcessTime = " << std::fixed << std::setprecision(10) << elapsed.count();
      std::cout << "\n";
    }
  }
};

}//end namespace
#endif



// private:

//   void multiRunSamplingForcingPeriod()
//   {
//     /* here we sample the forcing period,
//        which means that:
//        1. the material does not change
//        2. the other properties (like location) of the source do not change
//     */

//     /*
//      * create and store material prop
//      * only do it once since material does not change
//     */

//     auto matObj = createMaterialModel<sc_t>(parser_);
//     // compute rom Jacobians: only onces, these don't change
//     fom_t fomObj(meshInfo_);
//     fomObj.computeJacobiansWithMatProp(*matObj);
//     computeRomOperatorsUsingFomJacsWithMatProp(parser_, fomObj, phiVp_d_, phiSp_d_, Jvp_d_, Jsp_d_);

//     // create vector of signals
//     const auto periods = parser_.getValues(0);
//     using signal_t = Signal<sc_t>;
//     std::vector<signal_t> signals;
//     for (auto i=0; i<periods.size(); ++i)
//     {
//       signals.emplace_back(parser_.getSignal());
//       signals.back().resetPeriod(periods[i]);
//     }

//     std::cout << "Doing ROM with sampling of forcing period" << std::endl;
//     std::cout << "Total number of samples " << signals.size() << std::endl;

//     // create a forcing object with mem allocation
//     // (in loop below, only thing that changes is
//     // the signal NOT the location of the signal, so it is fine to
//     // create the nominal forcing and the in the loop below replace signal)
//     forcing_t forcingObj(parser_, meshInfo_, fomObj);

//     // index of where the source acts
//     const auto fIndex   = forcingObj.getVpGid();

//     // create vector containg phiVp^T*rhoInv which pre-multiplies the forcing
//     state_d_t phiVpRhoInvVec_d("phiVpRhoInv_d", nVp_);
//     // extract the target row from phiVp
//     auto phiVpRow_d = Kokkos::subview(phiVp_d_, fIndex, Kokkos::ALL());

//     // copy phiVp(fIndex,:) to phiVpRhoInvVec
//     Kokkos::deep_copy(phiVpRhoInvVec_d, phiVpRow_d);

//     // extract the target value from rhoInv
//     const auto rhoInvVpValue = fomObj.viewInvDensityHost(dofId::vp)(fIndex);
//     // scale
//     KokkosBlas::scal(phiVpRhoInvVec_d, rhoInvVpValue, phiVpRhoInvVec_d);

//     // loop over signals
//     for (std::size_t iSample=0; iSample<signals.size(); ++iSample)
//     {
//       // replace signal (no new allocations happen here)
//       forcingObj.replaceSignal(signals[iSample]);

//       // need to recheck that the new signal still meets conditions
//       doChecks(forcingObj, fomObj);

//       // reset observer
//       observerObj_.prepForNewRun(iSample);

//       // run
//       kokkosapp::runRomRankOneForcing(parser_.getNumSteps(), parser_.getTimeStepSize(),
// 				      phiVpRhoInvVec_d, forcingObj, observerObj_,
// 				      Jvp_d_, Jsp_d_, xVp_d_, xSp_d_);
//       processCollectedData(iSample);
//     }
//   }

//   template <typename forcing_t, typename fom_t>
//   void doChecks(const forcing_t & forcing, fom_t & fomObj){
//     if (parser_.checkDispersion())
//       checkDispersionCriterion(meshInfo_, forcing.getMaxFreq(),
//     			       fomObj.getMinShearWaveVelocity());

//     if (parser_.checkCfl())
//       checkCfl(meshInfo_, parser_.getTimeStepSize(), fomObj.getMaxShearWaveVelocity());
//   }

//   void processCollectedData(std::size_t iSample = 0)
//   {
//     // only process final/stored data if dummy basis = false
//     if (parser_.enableRandomDummyBasis() == false)
//     {

//       const auto startTime  = std::chrono::high_resolution_clock::now();

//       if(parser_.enableSnapshotMatrix()){
// 	observerObj_.writeSnapshotMatrixToFile(dofId::vp);
// 	observerObj_.writeSnapshotMatrixToFile(dofId::sp);
//       }

//       const auto finishTime = std::chrono::high_resolution_clock::now();
//       const std::chrono::duration<double> elapsed = finishTime - startTime;
//       std::cout << "\nfinalProcessTime = " << std::fixed << std::setprecision(10)
// 		<< elapsed.count();
//       std::cout << "\n";

//       // if(parser_.enableWriteFinalState()){
//       // 	processFinalStateToFile(dofId::vp, iSample);
//       // 	processFinalStateToFile(dofId::sp, iSample);
//       // }
//     }
//   }

//   // void processFinalStateToFile(const dofId dof, std::size_t iSample = 0)
//   // {
//   //   constexpr auto one	= constants<sc_t>::one();
//   //   constexpr auto zero	= constants<sc_t>::zero();

//   //   const auto dofName	  = dofIdToString(dof);
//   //   const auto romState_d = (dof==dofId::vp) ? xVp_d_ : xSp_d_;
//   //   const auto phi_d	  = (dof==dofId::vp) ? phiVp_d_ : phiSp_d_;

//   //   const auto fomStateSize     = (dof==dofId::vp) ? nVpFom_ : nSpFom_;
//   //   const auto fomStateFileName = "finalFomState_"+dofName+"_"+std::to_string(iSample);
//   //   const auto romStateFileName = "finalRomState_"+dofName+"_"+std::to_string(iSample);

//   //   state_d_t fomState_d("fomState_d", fomStateSize);
//   //   const char ct_N	= 'N';
//   //   KokkosBlas::gemv(&ct_N, one, phi_d, romState_d, zero, fomState_d);

//   //   auto fomState_h = Kokkos::create_mirror_view(fomState_d);
//   //   auto romState_h = Kokkos::create_mirror_view(romState_d);
//   //   Kokkos::deep_copy(fomState_h, fomState_d);
//   //   Kokkos::deep_copy(romState_h, romState_d);

//   //   writeToFile(fomStateFileName, fomState_h, parser_.writeFinalStateBinary());
//   //   // the rom state is small, so use ascii
//   //   writeToFile(romStateFileName, romState_h, false);
//   // }
