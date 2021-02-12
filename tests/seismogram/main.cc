
#include "./shared/all.hpp"
#include "./kokkos/types.hpp"
#include "./kokkos/shwavepp.hpp"

int main(int argc, char *argv[])
{
  Kokkos::initialize (argc, argv);
  {
    std::string sentinel = "PASS";

    using sc_t = typename kokkosapp::commonTypes::scalar_type;
    using parser_t = typename kokkosapp::commonTypes::parser_type;
    using mesh_info_t = typename kokkosapp::commonTypes::mesh_info_type;
    using seismo_t =  typename kokkosapp::commonTypes::seismogram_type;

    parser_t parser(argc, argv);
    mesh_info_t meshInfo(parser.getMeshDir());

    auto matObj = createMaterialModel<sc_t>(parser, meshInfo);
    kokkosapp::ShWavePP<kokkosapp::commonTypes> appObj(meshInfo, *matObj);

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
