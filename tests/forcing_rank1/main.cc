
#include "./fom/fom_problem_rank_one_forcing.hpp"

int main(int argc, char *argv[])
{
  Kokkos::initialize (argc, argv);
  {
    using sc_t = double;
    using parser_t = kokkosapp::commonTypes::parser_t;
    using prob_t = kokkosapp::FomProblemRankOneForcing;
    using app_t  = typename prob_t::fom_t;
    using mesh_info_t = typename prob_t::mesh_info_t;
    using forcing_t = typename prob_t::forcing_t;

    parser_t parser(argc, argv);
    mesh_info_t meshInfo(parser.getMeshDir());
    app_t appObj(meshInfo);

    auto matObj = createMaterialModel<sc_t>(parser);
    appObj.computeJacobians(*matObj);

    forcing_t forcingObj(parser, meshInfo, appObj);

    if (forcingObj.getMaxFreq() != 0.04){ std::puts("FAIL"); return 0; }
    if (forcingObj.getVpGid() != 636){ std::puts("FAIL"); return 0; }

    std::puts("PASS");
  }
  Kokkos::finalize();

  return 0;
}
