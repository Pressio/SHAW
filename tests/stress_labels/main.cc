
#include "./fom/fom_problem_rank_one_forcing.hpp"

int main(int argc, char *argv[])
{
  Kokkos::initialize (argc, argv);
  {
    std::string sentinel = "PASS";

    using sc_t = double;
    using parser_t = kokkosapp::commonTypes::parser_t;
    using prob_t = kokkosapp::FomProblemRankOneForcing;
    using app_t  = typename prob_t::fom_t;
    using mesh_info_t = typename prob_t::mesh_info_t;
    using forcing_t = typename prob_t::forcing_t;

    parser_t parser(argc, argv);
    mesh_info_t meshInfo(parser.getMeshDir());
    app_t appObj(meshInfo);

    auto matObj = createMaterialModel<sc_t>(parser, meshInfo);
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
