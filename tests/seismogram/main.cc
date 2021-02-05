
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
    using seismo_t = typename prob_t::seismogram_t;
    using mesh_info_t = typename prob_t::mesh_info_t;

    parser_t parser(argc, argv);
    mesh_info_t meshInfo(parser.getMeshDir());
    app_t appObj(meshInfo);

    auto matObj = createMaterialModel<sc_t>(parser, meshInfo);
    appObj.computeJacobians(*matObj);

    seismo_t seismoObj(parser, meshInfo, appObj);

    const auto gids = seismoObj.viewMappedGids();
    if (gids.extent(0) != 4){ sentinel = "FAIL"; }
    if (gids(0) != 1027){ sentinel = "FAIL"; }
    if (gids(1) != 1034){ sentinel = "FAIL"; }
    if (gids(2) != 1053){ sentinel = "FAIL"; }
    if (gids(3) != 1064){ sentinel = "FAIL"; }

    std::puts(sentinel.c_str());
  }
  Kokkos::finalize();

  return 0;
}
