
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

    if (meshInfo.getMeshDir() != "fullMesh21x51") sentinel = "FAIL";

    {
    const auto d1 = std::abs(meshInfo.getAngularSpacing() - 0.06283185307179586787);
    const auto d2 = std::abs(meshInfo.getRadialSpacing()  - 144550.00000000001136868);
    const auto d3 = std::abs(meshInfo.getAngularSpacingInverse() - 15.9154943092);
    const auto d4 = std::abs(meshInfo.getRadialSpacingInverse() - 0.000006918021446);
    if(d1>1e-10) sentinel="FAIL";
    if(d2>1e-10) sentinel="FAIL";
    if(d3>1e-10) sentinel="FAIL";
    if(d4>1e-10) sentinel="FAIL";
    }

    {
    const auto d1 = std::abs(meshInfo.getMinRadius() - 3480000.);
    const auto d2 = std::abs(meshInfo.getMinRadiusKm()  - 3480.);
    const auto d3 = std::abs(meshInfo.getMaxRadius() - 6371000.);
    const auto d4 = std::abs(meshInfo.getMaxRadiusKm()  - 6371.);
    if(d1>1e-10) sentinel="FAIL";
    if(d2>1e-10) sentinel="FAIL";
    if(d3>1e-10) sentinel="FAIL";
    if(d4>1e-10) sentinel="FAIL";
    }

    {
      const auto & db = meshInfo.viewDomainBounds();
      const auto d1 = std::abs(db[0] - 0.);
      const auto d2 = std::abs(db[1] - 180.);
      const auto d3 = std::abs(db[2] - 3480000.);
      const auto d4 = std::abs(db[3] - 6371000.);
      if(d1>1e-10) sentinel="FAIL";
      if(d2>1e-10) sentinel="FAIL";
      if(d3>1e-10) sentinel="FAIL";
      if(d4>1e-10) sentinel="FAIL";
    }

    {
      const auto d1 = std::abs(meshInfo.getMinArc() - 218654.8486898496);
      const auto d2 = std::abs(meshInfo.getMaxArc() - 400301.7359204115);
      if(d1>1e-10) sentinel="FAIL";
      if(d2>1e-10) sentinel="FAIL";
    }

    if (meshInfo.getNumPtsAlongR() != 21) sentinel = "FAIL";
    if (meshInfo.getNumPtsAlongTheta() != 51) sentinel = "FAIL";
    if (meshInfo.getNumVpPts() != 1071) sentinel = "FAIL";
    if (meshInfo.getNumSpPts() != 2070) sentinel = "FAIL";

    std::puts(sentinel.c_str());
  }
  Kokkos::finalize();

  return 0;
}
