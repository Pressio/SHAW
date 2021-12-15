
#include "CLI11.hpp"
#include "Kokkos_Core.hpp"
#include "../shared/constants.hpp"
#include "../shared/meta_kokkos.hpp"
#include "../shared/io/read_basis.hpp"
#include "../shared/io/matrix_write.hpp"
#include "../shared/io/matrix_read.hpp"
#include "../shared/io/vector_write.hpp"
#include "KokkosBlas2_gemv.hpp"
#include "utility"

using scalar_type = double;
using exe_space = Kokkos::DefaultHostExecutionSpace;
using kll = Kokkos::LayoutLeft;

struct CmdArgs
{
  using pair_t = std::pair<std::string, std::string>;

  pair_t podModes = {};
  pair_t romSnaps = {};
  std::size_t romSize = {};
  std::size_t fSize = {};
  std::vector<std::size_t> timeSteps = {};
  std::size_t samplingFreq = {};
  std::string outputFormat = {};
  std::string outFileAppend = {};

  std::string podFilePath;
  bool podIsBinary;
  std::string romSnapFile;
  bool romSnapBinary;

  // on input, we have the time steps so we need to convert
  // from time steps to indices of the snapsshopt matrix
  std::vector<int> targetIndices = {};

  CLI::App app;

  CmdArgs() = delete;
  CmdArgs(const std::string & name,
	  int argc, char *argv[])
    : app(name)
  {
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

    try{
      app.parse(argc, argv);
    }
    catch (...){

    }

    podFilePath   = std::get<0>(podModes);
    podIsBinary   = std::get<1>(podModes)=="binary";
    romSnapFile   = std::get<0>(romSnaps);
    romSnapBinary = std::get<1>(romSnaps)=="binary";
    for (auto it : timeSteps){
      targetIndices.push_back( it / samplingFreq );
    }
  }

  void describe() const{
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
  }

};

auto load_pod_modes(const CmdArgs & args)
{
  using pod_t = Kokkos::View<scalar_type**, kll, exe_space>;

  // create phi with 1 row because readBasis will resize rows
  pod_t phi("phi", 1, args.romSize);
  readBasis(args.podFilePath, phi, args.podIsBinary);
  std::cout << "POD modes size = "
	    << phi.extent(0) << " "
	    << phi.extent(1) << std::endl;

  double min ={};
  double max = {};
  for (std::size_t i=0; i<phi.extent(0); ++i){
    for (std::size_t j=0; j<phi.extent(1); ++j){
      min = std::min(min, phi(i,j));
      max = std::max(max, phi(i,j));
    }
  }
  std::cout << "POD modes MIN/MAX = "
	    << min << " " << max << std::endl;

  return phi;
}

template<class PhiType>
void execute_rank1(const CmdArgs & args, PhiType phi)
{
  // load ROM snaps
  using snap_t = Kokkos::View<scalar_type**, kll, exe_space>;
  snap_t snapsRom("snapsRom", 1, 1);
  if (args.romSnapBinary){
    fillMatrixFromBinary(args.romSnapFile, snapsRom, true);
  }
  else{
    fillMatrixFromAscii(args.romSnapFile, snapsRom, true);
  }
  std::cout << "ROM snap size: "
	    << snapsRom.extent(0) << " "
	    << snapsRom.extent(1) << std::endl;

  // *** reconstruct fom ***
  using state_t = Kokkos::View<scalar_type*, exe_space>;
  state_t fomState("fomState", phi.extent(0));
  const char ct_N	= 'N';

  for (std::size_t i=0; i<args.timeSteps.size(); ++i)
    {
      const auto thisTimeStep = args.timeSteps[i];
      auto romState = Kokkos::subview(snapsRom, Kokkos::ALL(), args.targetIndices[i]-1);
      writeToFile("ROMSTATE.txt", romState, false, true);

      KokkosBlas::gemv(&ct_N, 1, phi, romState, 0, fomState);
      auto fomStateFile = "fomReconstructedState_timestep_"+std::to_string(thisTimeStep);
      if (args.outFileAppend.empty() == false){
	fomStateFile += "_" + args.outFileAppend;
      }

      writeToFile(fomStateFile, fomState, (args.outputFormat=="binary"), true);
    }
}

template<class PhiType>
void execute_rank2(const CmdArgs & args, PhiType phi)
{
  throw std::runtime_error("reconstruction for rank2 to renable");

  // ////////////////////////////
  // // this is rank-2 case
  // ////////////////////////////

  // // *** load ROM states ***
  // using snap_t = Kokkos::View<scalar_type***, kll, exe_space>;
  // snap_t snapsRom("snapsRom", 1, 1, 1);
  // if (romSnapBinary) fillMatrixFromBinary(romSnapFile, snapsRom);
  // else
  // 	throw std::runtime_error("Rank-2 ROM snaps ascii not supported yet");
  // std::cout << "ROM snap size: "
  // 		<< snapsRom.extent(0) << " "
  // 		<< snapsRom.extent(1) << " "
  // 		<< snapsRom.extent(2) << std::endl;

  // // *** reconstruct fom state ***
  // using state_t = Kokkos::View<scalar_type*, exe_space>;
  // state_t fomState("fomState", phi.extent(0));
  // const char ct_N	= 'N';

  // for (std::size_t fId=0; fId<snapsRom.extent(2); ++fId)
  // {
  // 	const std::string fString = "f_"+std::to_string(fId);

  // 	for (std::size_t i=0; i<timeSteps.size(); ++i)
  // 	{
  // 	  const auto thisTimeStep = timeSteps[i];
  // 	  const std::string tsString = "timestep_"+std::to_string(thisTimeStep);
  // 	  auto romState = Kokkos::subview(snapsRom, Kokkos::ALL(), targetIndices[i]-1, fId);
  // 	  KokkosBlas::gemv(&ct_N, 1, phi, romState, 0, fomState);
  // 	  auto fomStateFile = "fomReconstructedState_"+tsString+"_"+fString;
  // 	  if (outFileAppend.empty() == false) fomStateFile += "_" + outFileAppend;

  // 	  writeToFile(fomStateFile, fomState, (outputFormat=="binary"), true);
  // 	}
  // }
}

int main(int argc, char *argv[])
{
  CmdArgs args("FOM state reconstructor", argc, argv);
  args.describe();

  Kokkos::initialize (argc, argv);
  {
    const auto phi = load_pod_modes(args);

    if (args.fSize == 1){
      execute_rank1(args, phi);
    }
    else if (args.fSize >= 2){
      execute_rank2(args, phi);
    }
  }
  Kokkos::finalize();

  std::cout << "success" << std::endl;
  return 0;
}
