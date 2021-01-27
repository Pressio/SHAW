
#include "CLI11.hpp"
#include "Kokkos_Core.hpp"
#include "../shared/constants.hpp"
#include "../shared/meta.hpp"
#include "../shared/io.hpp"
#include "KokkosBlas2_gemv.hpp"
#include "utility"

int main(int argc, char *argv[])
{
  CLI::App app{"FOM state reconstructor"};

  using pair_t = std::pair<std::string, std::string>;
  pair_t podModes = {};
  pair_t romSnaps = {};
  std::size_t romSize = {};
  std::size_t fSize = {};
  std::vector<std::size_t> timeSteps = {};
  std::size_t samplingFreq = {};
  std::string outputFormat = {};
  std::string outFileAppend = {};

  app.add_option("--romsize", romSize,
		 "ROM size")->required();

  app.add_option("--podmodes", podModes,
		 "Pair: fullpath_POD_modes binary/ascii")->required();

  app.add_option("--romsnaps", romSnaps,
		 "Pair: fullpath_ROM_snaps binary/ascii")->required();

  app.add_option("--samplingfreq", samplingFreq,
		 "Sampling freq used to save snapshot")->required();

  app.add_option("--fsize", fSize,
		 "Forcing size used to generate ROM snapshots")->required();

  app.add_option("--timesteps", timeSteps,
		 "Target time steps")->required();

  app.add_option("--outformat", outputFormat,
		 "Outputformat: binary/ascii")->required();

  app.add_option("--outfileappend", outFileAppend,
		 "String to append to output file");

  CLI11_PARSE(app, argc, argv);

  std::cout << "POD modes = "
	    << std::get<0>(podModes) << " "
	    << std::get<1>(podModes) << std::endl;
  std::cout << "ROM snaps = "
	    << std::get<0>(romSnaps) << " "
	    << std::get<1>(romSnaps) << std::endl;

  std::cout << "ROM size = " << romSize << std::endl;
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
  const auto podFilePath = std::get<0>(podModes);
  const bool podIsBinary = std::get<1>(podModes)=="binary";
  const auto romSnapFile   = std::get<0>(romSnaps);
  const bool romSnapBinary = std::get<1>(romSnaps)=="binary";

  // on input, we have the time steps so we need to convert
  // from time steps to indices of the snapsshopt matrix
  std::vector<int> targetIndices = {};
  for (auto it : timeSteps){
    targetIndices.push_back( it / samplingFreq );
  }

  //------------------------------------------------------------
  Kokkos::initialize (argc, argv);
  {
    using sc_t = double;
    using exe_space = Kokkos::DefaultExecutionSpace;
    using kll = Kokkos::LayoutLeft;
    static_assert(is_host_space<exe_space>::value, "");

    // *** load pod modes ***
    using pod_t = Kokkos::View<sc_t**, kll, exe_space>;
    auto phi = readBasis<sc_t,std::size_t,pod_t>(podFilePath, romSize, podIsBinary);
    std::cout << "POD modes size = "
	      << phi.extent(0) << " "
	      << phi.extent(1) << std::endl;

    if (fSize == 1)
    {
      ////////////////////////////
      // this is rank-1 case
      ////////////////////////////

      // *** load ROM states ***
      using snap_t = Kokkos::View<sc_t**, kll, exe_space>;
      snap_t snapsRom("snapsRom", 1, 1);
      if (romSnapBinary) readBinaryMatrixWithSize(romSnapFile, snapsRom);
      else readAsciiMatrixWithSize(romSnapFile, snapsRom);
      std::cout << "ROM snap size: "
		<< snapsRom.extent(0) << " "
		<< snapsRom.extent(1) << std::endl;

      // *** reconstruct fom ***
      using state_t = Kokkos::View<sc_t*, exe_space>;
      state_t fomState("fomState", phi.extent(0));
      const char ct_N	= 'N';

      for (std::size_t i=0; i<timeSteps.size(); ++i)
      {
	const auto thisTimeStep = timeSteps[i];
	auto romState = Kokkos::subview(snapsRom, Kokkos::ALL(), targetIndices[i]-1);
	KokkosBlas::gemv(&ct_N, 1, phi, romState, 0, fomState);
	auto fomStateFile = "fomReconstructedState_timestep_"+std::to_string(thisTimeStep);
	if (outFileAppend.empty() == false) fomStateFile += "_" + outFileAppend;

	writeToFile(fomStateFile, fomState, (outputFormat=="binary"), true);
      }
    }
    else if (fSize >= 2)
    {
      ////////////////////////////
      // this is rank-2 case
      ////////////////////////////

      // *** load ROM states ***
      using snap_t = Kokkos::View<sc_t***, kll, exe_space>;
      snap_t snapsRom("snapsRom", 1, 1, 1);
      if (romSnapBinary) readBinaryMatrixWithSize(romSnapFile, snapsRom);
      else
	throw std::runtime_error("Rank-2 ROM snaps ascii not supported yet");
      std::cout << "ROM snap size: "
		<< snapsRom.extent(0) << " "
		<< snapsRom.extent(1) << " "
		<< snapsRom.extent(2) << std::endl;

      // *** reconstruct fom state ***
      using state_t = Kokkos::View<sc_t*, exe_space>;
      state_t fomState("fomState", phi.extent(0));
      const char ct_N	= 'N';

      for (std::size_t fId=0; fId<snapsRom.extent(2); ++fId)
      {
	const std::string fString = "f_"+std::to_string(fId);

	for (std::size_t i=0; i<timeSteps.size(); ++i)
	{
	  const auto thisTimeStep = timeSteps[i];
	  const std::string tsString = "timestep_"+std::to_string(thisTimeStep);
      	  auto romState = Kokkos::subview(snapsRom, Kokkos::ALL(), targetIndices[i]-1, fId);
      	  KokkosBlas::gemv(&ct_N, 1, phi, romState, 0, fomState);
      	  auto fomStateFile = "fomReconstructedState_"+tsString+"_"+fString;
	  if (outFileAppend.empty() == false) fomStateFile += "_" + outFileAppend;

      	  writeToFile(fomStateFile, fomState, (outputFormat=="binary"), true);
      	}
      }
    }
  }
  Kokkos::finalize();
  std::cout << "success" << std::endl;

  return 0;
}
