
#include "Eigen/Dense"
#include "Kokkos_Core.hpp"
#include <Kokkos_Random.hpp>
#include <KokkosBlas1_axpby.hpp>
#include <KokkosBlas1_fill.hpp>
#include "KokkosBlas2_gemv.hpp"
#include "../shared/constants.hpp"
#include "../shared/complexity.hpp"
#include "../shared/print_perf.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <vector>

namespace
{
using sc_t	    = double;

using klr	  = Kokkos::LayoutRight;
using kll	  = Kokkos::LayoutLeft;
using exe_space	  = Kokkos::DefaultExecutionSpace;
using ko_mat_ll_t = Kokkos::View<sc_t**, kll, exe_space>;
using ko_mat_lr_t = Kokkos::View<sc_t**, klr, exe_space>;
using ko_vec_t	  = Kokkos::View<sc_t*, exe_space>;
}//end anonym namespace

void gemvKokkos(int m, int p, int rep)
{
  // y = beta*y + alpha*A*x

  const char ct_N	= 'N';
  constexpr auto zero	= constants<sc_t>::zero();
  constexpr auto one	= constants<sc_t>::one();
  constexpr auto two	= constants<sc_t>::two();

  Kokkos::Random_XorShift64_Pool<exe_space> rand_pool(13718);
  ko_mat_lr_t A("A", m, p);
  ko_vec_t y("y", m);
  ko_vec_t x("x", p);
  Kokkos::fill_random(A, rand_pool, sc_t(10));
  Kokkos::fill_random(x, rand_pool, sc_t(10));
  KokkosBlas::fill(y, one);

  // warmup
  //KokkosBlas::gemm(&ct_N, &ct_N, one, A, B, one, C);

  std::array<double, 3> perfTimes = {1e32,0.,0.}; //min, max, ave
  const auto startTime = std::chrono::high_resolution_clock::now();
  for (auto i=0; i<rep; ++i){
    Kokkos::Timer timer;

    KokkosBlas::gemv(&ct_N, two, A, x, one, y);

    const double time = timer.seconds();
    perfTimes[0] = std::min(perfTimes[0], time);
    perfTimes[1] = std::max(perfTimes[1], time);
    perfTimes[2] += time;
  }
  const auto finishTime = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finishTime - startTime;
  std::cout << "Elapsed time: "
	    << std::fixed << std::setprecision(10)
	    << elapsed.count() << std::endl;

  double memCostMB, flopsCost = 0.;
  Complexity<sc_t>::gemv_beta_one(m, p, memCostMB, flopsCost);
  kokkosapp::printPerf(rep, perfTimes, memCostMB, flopsCost);
}

int main(int argc, char *argv[])
{
  // *** parse args ***
  int m = 0;
  int p = 0;
  int rep  = 0;

  if(argc>=2){
    m = std::atoi(argv[1]);
    p = std::atoi(argv[2]);
    rep  = std::atoi(argv[3]);
  }

  std::cout << m << " "
	    << p << " "
	    << rep << std::endl;

  Kokkos::initialize (argc, argv);
  gemvKokkos(m, p, rep);
  Kokkos::finalize();

  return 0;
}
