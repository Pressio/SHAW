/*
//@HEADER
// ************************************************************************
//
// seismogram.hpp
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

#ifndef SHAXIPP_SEISMOGRAM_HPP_
#define SHAXIPP_SEISMOGRAM_HPP_

template <typename gids_t, typename state_t, typename dest_t>
struct CopySeis
{
  std::size_t count_;
  gids_t gids_;
  state_t x_;
  dest_t M_;

  CopySeis(std::size_t count, gids_t gids, state_t x, dest_t M)
    : count_(count), gids_(gids), x_(x), M_(M){}

  template <typename _state_t = state_t>
  KOKKOS_INLINE_FUNCTION
  typename std::enable_if<is_kokkos_1dview<_state_t>::value>::type
  operator() (const std::size_t & i) const
  {
    M_(i, count_, 0) = x_(gids_(i));
  }

  template <typename _state_t = state_t>
  KOKKOS_INLINE_FUNCTION
  typename std::enable_if<is_kokkos_2dview<_state_t>::value>::type
  operator() (const std::size_t & i) const
  {
    for (std::size_t j=0; j<M_.extent(2); ++j)
      M_(i, count_, j) = x_(gids_(i), j);
  }
};


template <typename scalar_t>
class Seismogram
{
  // container type to store the data
  // here we have to specify the layout because of how we write to file
  using matrix_t = Kokkos::View<scalar_t***, Kokkos::LayoutLeft, Kokkos::HostSpace>;

  // type to store global ids
  using gids_t = Kokkos::View<std::size_t*, Kokkos::HostSpace>;

  // flag to enable/disable collection of seismogram data
  bool enable_ = {};

  // whether we need to store in binary
  bool useBinaryIO_   = {};

  // filename to store seigmogram data
  std::string seismoFileName_ = {};

  // sampling frequency
  std::size_t freq_  = 0;

  std::size_t count_ = 0;

  // number of receivers
  std::size_t numReceivers_ {};

  // list of velocity gids identifying the elements in the velocity state
  // where we need to sample at
  gids_t targetGids_;

  // data matrix: each row identifies a receiver, the columns store vp(t)
  matrix_t MM_;

  // runID used when we run many samples to prepend file
  std::size_t runID_ = 0;

public:
  template <typename parser_t, typename mesh_info_t, typename app_t>
  Seismogram(const parser_t & parser,
	     const mesh_info_t & meshInfo,
	     const app_t & appObj,
	     std::size_t fSize = 1)
    : enable_{parser.enableSeismogram()},
      useBinaryIO_(parser.writeSeismogramBinary()),
      seismoFileName_{parser.getSeismogramFileName()}
  {
    if (enable_){
      std::cout << std::endl;
      std::cout << "*** Constructing seismogram ***" << std::endl;

      // all receivers are located on surface
      const auto receiverRadM = meshInfo.viewDomainBounds()[3];

      // get list of nominal angles of the receivers
      const auto nominalAngles = parser.getSeismoReceiversAnglesDeg();
      numReceivers_ = nominalAngles.size();

      // loop over nominal locations and map them to grid
      // this is because nominal values are typicall round angles like 30, 60,
      // but the grid does not have any point exactly at those locations
      // so we need to map the nominal angles to grid points that are closest
      // to the desired values
      Kokkos::resize(targetGids_, numReceivers_);
      for (auto i=0; i<numReceivers_; ++i)
      {
	const auto thisAngle = nominalAngles[i];
	targetGids_(i) = mapNominalLocationToVelocityGridPoint
	  (thisAngle,
	   receiverRadM,
	   meshInfo.viewDomainBounds(),
	   meshInfo.getNumVpPts(),
	   appObj.viewGidListHost(dofId::vp),
	   appObj.viewCoordsHost(dofId::vp));
      }
      std::cout << "Done Mapping receivers to grid " << std::endl;

      freq_ = parser.getSeismoFreq();
      const auto Nsteps = parser.getNumSteps();
      // make sure number of steps is divisible by sampling frequency
      if ( Nsteps % freq_ == 0){
	const auto numCols = Nsteps/freq_;
	Kokkos::resize(MM_, numReceivers_, numCols, fSize);
      }
      else{
	throw std::runtime_error("Seismogram sampling frequency not a divisor of steps");
      }

      // estimate how much memory is needed to store seismogram
      const double mem = MM_.extent(0)*MM_.extent(1)*MM_.extent(2) * sizeof(scalar_t);
      std::cout << "Seismogram [GB] = " << mem/(1024.*1024.*1024.) << std::endl;
      std::cout << std::endl;
    }
  }

  bool enabled() const
  {
    return enable_;
  }

  const gids_t & viewMappedGids() const{
    return targetGids_;
  }

  void prepForNewRun(const std::size_t & sampleID){
    // assumes the new run has same sampling frequncies as before
    count_ = {0};
    runID_ = sampleID;
  }

  template <typename state_t>
  void storeVelocitySignalAtReceivers(std::size_t step,
				      const state_t & xhv)
  {
    static_assert
      (std::is_same<typename state_t::memory_space, Kokkos::HostSpace>::value,
       "View not accessible on host");

    if (enable_)
    {
      if ( step % freq_ == 0 and step > 0){
	using functor_t = CopySeis<gids_t, state_t, matrix_t>;
	functor_t fnc(count_, targetGids_, xhv, MM_);

	// must specify an host exespace here otherwise it picks the default
	// which might be a device one
	using copy_exespace = Kokkos::DefaultHostExecutionSpace;
	Kokkos::RangePolicy<copy_exespace> policy(0, targetGids_.extent(0));
	Kokkos::parallel_for(policy, fnc);
	count_++;
      }
    }
  }

  void writeReceiversToFile() const
  {
    if (enable_){
      std::cout << "Writing seismogram ";
      auto fN2 = seismoFileName_ + "_" + std::to_string(runID_);
      if (MM_.extent(2) == 1){
	const auto Mv = Kokkos::subview(MM_, Kokkos::ALL(), Kokkos::ALL(), 0);
	writeToFile(fN2, Mv, useBinaryIO_, false); //false to not print size
      }
      else{
	writeToFile(fN2, MM_, useBinaryIO_, false); //false to not print size
      }
      std::cout << "... Done" << std::endl;
    }
  }
};

#endif
