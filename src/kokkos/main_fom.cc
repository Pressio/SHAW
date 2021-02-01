
#include "./fom/fom_problem_rank_one_forcing.hpp"
#include "./fom/fom_problem_rank_two_forcing.hpp"

int main(int argc, char *argv[])
{
  Kokkos::initialize (argc, argv);
  {
    using parser_t = kokkosapp::commonTypes::parser_t;
    parser_t parser(argc, argv);

    if(parser.enableForcingBatching())
    {
      using prob_t = kokkosapp::FomProblemRankTwoForcing;
      prob_t problem(parser);
      problem.execute();
    }
    else{
      using prob_t = kokkosapp::FomProblemRankOneForcing;
      prob_t problem(parser);
      problem.execute();
    }
  }
  Kokkos::finalize();

  return 0;
}
