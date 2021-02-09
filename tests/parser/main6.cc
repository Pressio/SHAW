
#include "./shared/all.hpp"
#include "./kokkos/common_types.hpp"

int main(int argc, char *argv[])
{
  Kokkos::initialize (argc, argv);
  {
    using parser_t = kokkosapp::commonTypes::parser_type;
    parser_t parser(argc, argv);

    std::vector<bool> vb;

    // forcing
    vb.push_back(parser.getSourceSignalKind() == signalKind::ricker);
    vb.push_back(parser.viewAngles()[0]== 88.);
    vb.push_back(parser.viewDepths()[0] == 1122.);
    vb.push_back(parser.viewDelays()[0]== 10.);

    const auto & vals = parser.viewPeriods();
    vb.push_back(vals[0] == 40.);
    vb.push_back(vals[1] == 55.);
    vb.push_back(vals[2] == 65.);

    if (std::none_of(vb.begin(), vb.end(), std::logical_not<bool>()))
      std::puts("PASS");
    else
      std::puts("FAIL");
  }
  Kokkos::finalize();

  return 0;
}
