
#ifndef SHAXIPP_OBSERVER_HPP_
#define SHAXIPP_OBSERVER_HPP_

template <typename state_t, typename dest_t>
struct Copy
{
  std::size_t colIndex_;
  state_t x_;
  dest_t M_;

  Copy(std::size_t colIndex, state_t x, dest_t M)
    : colIndex_(colIndex), x_(x), M_(M){}

  template <typename _state_t = state_t>
  KOKKOS_INLINE_FUNCTION
  typename std::enable_if<is_kokkos_1dview<_state_t>::value>::type
  operator() (const std::size_t & i) const
  {
    M_(i, colIndex_, 0) = x_(i);
  }

  template <typename _state_t = state_t>
  KOKKOS_INLINE_FUNCTION
  typename std::enable_if<is_kokkos_2dview<_state_t>::value>::type
  operator() (const std::size_t & i) const
  {
    for (std::size_t j=0; j<M_.extent(2); ++j)
      M_(i, colIndex_, j) = x_(i,j);
  }
};


template <typename int_t, typename scalar_t>
struct StateObserver
{
  using matrix_t = Kokkos::View<scalar_t***, Kokkos::LayoutLeft, Kokkos::HostSpace>;

private:
  bool useBinaryIO_   = {};
  bool enableSnapMat_ = {};
  std::array<std::string,2> snapFileName_ = {};

  std::array<int_t, 2> numDofs_ = {};

  // to count the snapshots
  std::array<int_t, 2> count_ = {};
  // sampling frequency
  std::array<int_t, 2> snapshotFreq_ = {};

  // snapshot matrices
  matrix_t Avp_;
  matrix_t Asp_;

  // runID used when we run many samples to prepend file
  int_t runID_ = 0;

public:
  template <typename parser_t>
  StateObserver(int_t numDof_vp,
		int_t numDof_sp,
		const parser_t & parser,
		int_t fSize = 1)
    : useBinaryIO_(parser.writeSnapshotsBinary()),
      enableSnapMat_{parser.enableSnapshotMatrix()},
      snapFileName_{{parser.getSnapshotFileName(dofId::vp),
		     parser.getSnapshotFileName(dofId::sp)}},
      numDofs_{{numDof_vp, numDof_sp}},
      snapshotFreq_{{parser.getSnapshotFreq(dofId::vp),
		     parser.getSnapshotFreq(dofId::sp)}}
  {
    if (enableSnapMat_){
      std::cout << "\n";
      std::cout << "*** Constructing observer ***" << std::endl;

      const auto Nsteps = parser.getNumSteps();
      int_t numColsVp = 0;
      int_t numColsSp = 0;

      // make sure number of steps is divisible by sampling frequency
      if ( Nsteps % snapshotFreq_[0] == 0){
	numColsVp = Nsteps/snapshotFreq_[0];
      }
      else{
	throw std::runtime_error("Snapshot Vp frequency not a divisor of steps");
      }

      // make sure number of steps is divisible by sampling frequency
      if ( Nsteps % snapshotFreq_[1] == 0 ){
	numColsSp = Nsteps/snapshotFreq_[1];
      }
      else{
	throw std::runtime_error("Snapshot Sp frequency not a divisor of steps");
      }

      //resize matrix
      Kokkos::resize(Avp_, numDof_vp, numColsVp, fSize);
      Kokkos::resize(Asp_, numDof_sp, numColsSp, fSize);

      const double memAvp = Avp_.extent(0)*Avp_.extent(1)*Avp_.extent(2) * sizeof(scalar_t);
      const double memAsp = Asp_.extent(0)*Asp_.extent(1)*Asp_.extent(2) * sizeof(scalar_t);
      std::cout << "Observer: Vp snaps [GB] = " << memAvp/(1024.*1024.*1024.) << std::endl;
      std::cout << "Observer: Sp snaps [GB] = " << memAsp/(1024.*1024.*1024.) << std::endl;
    }
  }

  void prepForNewRun(int_t runIdIn)
  {
    // assumes the new run has same sampling frequncies as before
    count_ = {0,0};
    runID_ = runIdIn;
  }

  const auto & viewSnapshotMatrix(const dofId dof) const
  {
    switch (dof){
    case dofId::vp: return Avp_;
    case dofId::sp: return Asp_;
    }
  }

  template<typename state_t>
  void observe(const dofId dof, int_t step, const state_t x)
  {

    if (enableSnapMat_)
    {
      auto & A	      = (dof==dofId::vp) ? Avp_ : Asp_;
      const auto freq = (dof==dofId::vp) ? snapshotFreq_[0] : snapshotFreq_[1];
      auto & count    = (dof==dofId::vp) ? count_[0] : count_[1];

      if ( step % freq == 0 and step > 0)
      {
  	auto xhv = Kokkos::create_mirror_view(x);
	Kokkos::deep_copy(xhv, x);

	// might want to make the copy a bit more efficient
	using state_h_t = typename state_t::HostMirror;
	Copy<state_h_t, matrix_t> fnc(count, xhv, A);
	Kokkos::parallel_for(xhv.extent(0), fnc);
  	count++;
      }
    }
  }

  void writeSnapshotMatrixToFile(const dofId dof) const
  {
    if (enableSnapMat_)
    {
      std::cout << "Writing snapshots " + dofIdToString(dof);

      const auto dofName = dofIdToString(dof);
      auto & A = (dof==dofId::vp) ? Avp_ : Asp_;
      auto & fN = (dof==dofId::vp) ? snapFileName_[0] : snapFileName_[1];
      // append the runID to the file
      auto fN2 = fN + "_" + std::to_string(runID_);

      if (A.extent(2) == 1)
      {
	// if the forcing has rank-1 the third dim of A is 1
	const auto Av = Kokkos::subview(A, Kokkos::ALL(), Kokkos::ALL(), 0);
	writeToFile(fN2, Av, useBinaryIO_);
      }
      else{
	writeToFile(fN2, A, useBinaryIO_);
      }
      std::cout << "... Done" << std::endl;
    }
  }

};
#endif

// std::string finalFileName = "final_state_" + dofName;
// if runID != -1 it means we are doing multiple runs, so modify file name
// 		    finalFileName += "_" + std::to_string(runID_);

// if (enableSnapViz_){
//   const auto vizFreq = (dof==dofId::vp) ? vizFreq_[0] : vizFreq_[1];
//   if ( step % vizFreq == 0 and step > 0){
// 	const auto dofName = dofIdToString(dof);
// 	auto xhv = Kokkos::create_mirror_view(x);
// 	Kokkos::deep_copy(xhv, x);
//   	auto fileName = "state_" + dofName + "_" + std::to_string(step);
//   	if (runID_ != -1) fileName += "_" + std::to_string(runID_);
//   	writeMatrixWithSizeToFile(fileName, xhv, useBinaryIO_);
//   }
// }
