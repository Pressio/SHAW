
#include "./shared/all.hpp"
#include "./kokkos/common_types.hpp"
#include "./kokkos/shwavepp.hpp"

template<class T>
std::string checkVpGraph(const T & appObj)
{
  std::string sentinel = "PASS";

  // check graph for Vp
  const std::vector<std::vector<std::size_t>> goldGVp =
    {
      {1, 0, 1, 0},
      {1, 2, 3, 2},
      {3, 4, 5, 4},
      {5, 6, 7, 6},
      {7, 8, 7, 8},
      {10,  9, 10,  0},
      {10, 11, 12,  2},
      {12, 13, 14,  4},
      {14, 15, 16,  6},
      {16, 17, 16,  8},
      {19, 18, 19,  9},
      {19, 20, 21, 11},
      {21, 22, 23, 13},
      {23, 24, 25, 15},
      {25, 26, 25, 17},
      {27, 18, 27, 18},
      {27, 20, 28, 20},
      {28, 22, 29, 22},
      {29, 24, 30, 24},
      {30, 26, 30, 26}
    };

  auto Gvp = appObj.viewVelocityGraphHost();
  if (Gvp.extent(0) != goldGVp.size()) sentinel = "FAIL";
  if (Gvp.extent(1) != 5) sentinel = "FAIL";
  for (std::size_t i=0; i<Gvp.extent(0); ++i)
    {
      if (Gvp(i,0) != i){
	sentinel = "FAIL";
	break;
      }

      for (std::size_t j=1; j<Gvp.extent(1); ++j){
	if (Gvp(i,j) != goldGVp[i][j-1]){
	  sentinel = "FAIL";
	  break;
	}
      }
    }

  return sentinel;
}

template<class T>
std::string checkSpGraph(const T & appObj)
{
  std::string sentinel = "PASS";

  // check graph for Sp
  const std::vector<std::vector<std::size_t>> goldGSp =
    {
      {5, 0},
      {0, 1},
      {6, 1},
      {1, 2},
      {7, 2},
      {2, 3},
      {8, 3},
      {3, 4},
      {9, 4},
      {10 , 5},
      {5, 6},
      {11 , 6},
      {6, 7},
      {12 , 7},
      {7, 8},
      {13 , 8},
      {8, 9},
      {14 , 9},
      {15, 10},
      {10, 11},
      {16, 11},
      {11, 12},
      {17, 12},
      {12, 13},
      {18, 13},
      {13, 14},
      {19, 14},
      {15, 16},
      {16, 17},
      {17, 18},
      {18, 19}
    };

  auto Gvp = appObj.viewStressGraphHost();
  if (Gvp.extent(0) != goldGSp.size()) sentinel = "FAIL";
  if (Gvp.extent(1) != 3) sentinel = "FAIL";
  for (std::size_t i=0; i<Gvp.extent(0); ++i)
    {
      if (Gvp(i,0) != i){
	sentinel = "FAIL";
	break;
      }

      for (std::size_t j=1; j<Gvp.extent(1); ++j){
	if (Gvp(i,j) != goldGSp[i][j-1]){
	  sentinel = "FAIL";
	  break;
	}
      }
    }

  return sentinel;
}


int main(int argc, char *argv[])
{
  Kokkos::initialize (argc, argv);
  {
    using sc_t = typename kokkosapp::commonTypes::scalar_type;
    using parser_t = typename kokkosapp::commonTypes::parser_type;
    using mesh_info_t = typename kokkosapp::commonTypes::mesh_info_type;

    parser_t parser(argc, argv);
    mesh_info_t meshInfo(parser.getMeshDir());
    auto matObj = createMaterialModel<sc_t>(parser, meshInfo);
    kokkosapp::ShWavePP<kokkosapp::commonTypes> appObj(meshInfo);
    appObj.computeJacobians(*matObj);

    const std::string sentinelVp = checkVpGraph(appObj);
    const std::string sentinelSp = checkSpGraph(appObj);

    if (sentinelVp == "PASS" and sentinelSp=="PASS"){
      std::puts("PASS");
    }
    else{
      std::puts("FAIL");
    }
  }
  Kokkos::finalize();

  return 0;
}
