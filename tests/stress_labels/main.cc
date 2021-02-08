
#include "./shared/all.hpp"
#include "./kokkos/common_types.hpp"
#include "./kokkos/shwavepp.hpp"

int main(int argc, char *argv[])
{
  Kokkos::initialize (argc, argv);
  {
    std::string sentinel = "PASS";

    using sc_t = typename kokkosapp::commonTypes::scalar_type;
    using parser_t = typename kokkosapp::commonTypes::parser_type;
    using mesh_info_t = typename kokkosapp::commonTypes::mesh_info_type;

    parser_t parser(argc, argv);
    mesh_info_t meshInfo(parser.getMeshDir());

    auto matObj = createMaterialModel<sc_t>(parser, meshInfo);
    kokkosapp::ShWavePP<kokkosapp::commonTypes> appObj(meshInfo);
    appObj.computeJacobians(*matObj);

    const auto labels = appObj.viewLabelsHost(dofId::sp);

    // gold vector with labels
    const std::vector<int> goldLabels =
      {1,2,1,2,1,2,1,2,1,1,2,1,2,1,2,1,2,1,1,2,1,2,1,2,1,2,1,2,2,2,2};

    // *** do check ***
    if (labels.extent(0) != goldLabels.size()){
      sentinel = "FAIL";
    }

    for (std::size_t i=0; i<labels.size(); ++i){
      if (labels(i) != goldLabels[i]){
	sentinel = "FAIL";
      }
    }

    std::puts(sentinel.c_str());
  }
  Kokkos::finalize();

  return 0;
}
