
#include "./fom/fom_problem_rank_one_forcing.hpp"
#include "./fom/fom_problem_rank_two_forcing.hpp"

template<typename parser_t>
std::shared_ptr<kokkosapp::FomProblemBase> createFomProblem(const parser_t & parser)
{
  if(parser.enableForcingBatching()){
    using ret_t = kokkosapp::FomProblemRankTwoForcing;
    return std::make_shared<ret_t>(parser);
  }
  else{
    using ret_t = kokkosapp::FomProblemRankOneForcing;
    return std::make_shared<ret_t>(parser);
  }
}

int main(int argc, char *argv[])
{
  Kokkos::initialize (argc, argv);
  {
    using parser_t = kokkosapp::commonTypes::parser_t;
    parser_t parser(argc, argv);

    auto problem = createFomProblem(parser);
    problem->execute();
  }
  Kokkos::finalize();

  return 0;
}
