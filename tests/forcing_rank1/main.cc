
#include "./shared/all.hpp"
#include "./kokkos/types.hpp"
#include "./kokkos/shwavepp.hpp"

int main(int argc, char *argv[])
{
  Kokkos::initialize (argc, argv);
  {
    using types = kokkosapp::rank1Types;

    using sc_t = typename types::scalar_type;
    using parser_t = typename types::parser_type;
    using mesh_info_t = typename types::mesh_info_type;
    using forcing_t = typename types::forcing_type;

    parser_t parser(argc, argv);
    mesh_info_t meshInfo(parser.getMeshDir());

    auto matObj = createMaterialModel<sc_t>(parser, meshInfo);
    kokkosapp::ShWavePP<types> appObj(meshInfo, *matObj);
    forcing_t forcingObj(parser, meshInfo, appObj);

    if (forcingObj.getMaxFreq() != 0.04){ std::puts("FAIL"); return 0; }
    if (forcingObj.getVpGid() != 636){ std::puts("FAIL"); return 0; }

    std::puts("PASS");
  }
  Kokkos::finalize();

  return 0;
}
