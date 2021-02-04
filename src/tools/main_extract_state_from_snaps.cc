
#include "CLI11.hpp"
#include "Kokkos_Core.hpp"
#include "KokkosBlas2_gemv.hpp"
#include "../shared/constants.hpp"
#include "../shared/meta/meta_eigen.hpp"
#include "../shared/meta/meta_kokkos.hpp"
#include "../shared/io/matrix_write.hpp"
#include "../shared/io/matrix_read.hpp"
#include "../shared/io/read_basis.hpp"
#include "../shared/io/read_reference_state.hpp"
#include "../shared/io/vector_write.hpp"
#include "../shared/io/vector_read.hpp"
#include "utility"

int main(int argc, char *argv[])
{
  CLI::App app{"Extract state(s) at target time steps from a snapshot matrix"};

  using pair_t = std::pair<std::string, std::string>;
  pair_t snaps = {};
  std::size_t fSize = 1;
  std::vector<std::size_t> timeSteps = {};
  std::size_t samplingFreq = {};
  std::string outputFormat = {};
  std::string outFileAppend = {};

  app.add_option("--snaps", snaps,
		 "Pair: fullpath_to_snaps binary/ascii")->required();

  app.add_option("--samplingfreq", samplingFreq,
		 "Sampling freq used to save snapshot")->required();

  app.add_option("--fsize", fSize,
		 "Forcing size used to generate snapshots")->required();

  app.add_option("--timesteps", timeSteps,
		 "Target time steps to extract")->required();

  app.add_option("--outformat", outputFormat,
		 "Outputformat: binary/ascii")->required();

  app.add_option("--outfileappend", outFileAppend,
		 "String to append to output file");

  CLI11_PARSE(app, argc, argv);

  std::cout << "snaps file = " << std::get<0>(snaps) << std::endl;
  std::cout << "snaps format = " << std::get<1>(snaps) << std::endl;
  std::cout << "Size of f = " << fSize << std::endl;

  if (timeSteps[0]==-1){
    std::cout << "Target time steps = only last step";
  }
  else{
    std::cout << "Target time steps = ";
    for (auto & it : timeSteps)
      std::cout << it << " ";
    std::cout << std::endl;
  }
  std::cout << "Sampling freq = " << samplingFreq << std::endl;
  std::cout << "Output format = " << outputFormat << std::endl;

  //------------------------------------------------------------
  const auto snapFile   = std::get<0>(snaps);
  const bool snapBinary = std::get<1>(snaps)=="binary";

  // on input, we have the time steps so we need to convert
  // from time steps to indices of the snapsshopt matrix
  std::vector<int> targetIndices = {};
  for (auto it : timeSteps){
    targetIndices.push_back( it / samplingFreq );
  }

  //------------------------------------------------------------
  Kokkos::initialize(); //argc, argv);
  {
    using sc_t = double;
    using exe_space = Kokkos::DefaultExecutionSpace;
    using kll = Kokkos::LayoutLeft;
    static_assert(is_host_space<exe_space>::value, "");

    if (fSize == 1)
    {
      ////////////////////////////
      // this is rank-1 case
      ////////////////////////////

      // *** load states ***
      using snap_t = Kokkos::View<sc_t**, kll, exe_space>;
      snap_t snaps("snaps", 1, 1);
      if (snapBinary) readBinaryMatrixWithSize(snapFile, snaps);
      else readAsciiMatrixWithSize(snapFile, snaps);
      std::cout << "snap size: "
		<< snaps.extent(0) << " "
		<< snaps.extent(1) << std::endl;

      for (std::size_t i=0; i<timeSteps.size(); ++i){
	const auto thisTimeStep = timeSteps[i];
	auto stateFile = "state_timestep_" + std::to_string(thisTimeStep);
	if (outFileAppend.empty() == false) stateFile += "_" + outFileAppend;

	auto state = Kokkos::subview(snaps, Kokkos::ALL(), targetIndices[i]-1);
	writeToFile(stateFile, state, (outputFormat=="binary"), true);
      }
    }
    else if (fSize >= 2)
    {
      ////////////////////////////
      // this is rank-2 case
      ////////////////////////////

      // *** load states ***
      using snap_t = Kokkos::View<sc_t***, kll, exe_space>;
      snap_t snaps("snaps", 1, 1, 1);
      if (snapBinary) readBinaryMatrixWithSize(snapFile, snaps);
      else
	throw std::runtime_error("Rank-2  snaps ascii not supported yet");
      std::cout << "snap size: "
		<< snaps.extent(0) << " "
		<< snaps.extent(1) << " "
		<< snaps.extent(2) << std::endl;

      for (std::size_t fId=0; fId<snaps.extent(2); ++fId)
      {
	const std::string fString = "f_"+std::to_string(fId);

	for (std::size_t i=0; i<timeSteps.size(); ++i)
	{
	  const auto thisTimeStep = timeSteps[i];
	  const auto tsString = "timestep_"+std::to_string(thisTimeStep);
      	  auto stateFile = "state_"+tsString+"_"+fString;
	  if (outFileAppend.empty() == false) stateFile += "_" + outFileAppend;

      	  auto state = Kokkos::subview(snaps, Kokkos::ALL(), targetIndices[i]-1, fId);
      	  writeToFile(stateFile, state, (outputFormat=="binary"), true);
      	}
      }
    }
  }
  Kokkos::finalize();

  std::cout << "success" << std::endl;

  return 0;
}
