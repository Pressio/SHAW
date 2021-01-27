
#include "Eigen/Dense"
#include "Kokkos_Core.hpp"
#include <Kokkos_Random.hpp>
#include <KokkosBlas1_axpby.hpp>
#include <KokkosBlas1_fill.hpp>
#include "KokkosBlas3_gemm.hpp"
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

using eig_mat_ll_t  = Eigen::Matrix<sc_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
using eig_mat_lr_t  = Eigen::Matrix<sc_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using eig_vec_t     = Eigen::Matrix<sc_t, Eigen::Dynamic, 1>;

using klr	  = Kokkos::LayoutRight;
using kll	  = Kokkos::LayoutLeft;
using exe_space	  = Kokkos::DefaultExecutionSpace;
using ko_mat_ll_t = Kokkos::View<sc_t**, kll, exe_space>;
using ko_mat_lr_t = Kokkos::View<sc_t**, klr, exe_space>;
}//end anonym namespace

void gemmEigen(int m, int p, int n, int rep)
{
  eig_mat_ll_t A(m, p); A.setRandom();
  eig_mat_ll_t B(p, n); B.setRandom();
  eig_mat_ll_t C(m, n); C.setOnes();

  // first touch, warmup
  C.noalias() += A * B;

  const auto startTime = std::chrono::high_resolution_clock::now();

  for (auto i=0; i<rep; ++i){
    C.noalias() += A * B;
  }

  const auto finishTime = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finishTime - startTime;
  std::cout << "Elapsed time: "
	    << std::fixed << std::setprecision(10)
	    << elapsed.count() << std::endl;

  const sc_t gFlopCount = double(m)*n*p*rep*2*1e-9;
  std::cout << gFlopCount/elapsed.count()  <<  " GFLOPS \n";
}

void gemmKokkos(int m, int p, int n, int rep)
{
  // C = beta*C + alpha*A*B

  const char ct_N	= 'N';
  constexpr auto zero	= constants<sc_t>::zero();
  constexpr auto one	= constants<sc_t>::one();

  Kokkos::Random_XorShift64_Pool<exe_space> rand_pool(13718);
  ko_mat_lr_t A("A", m, p);
  ko_mat_lr_t B("B", p, n);
  ko_mat_lr_t C("C", m, n);
  Kokkos::fill_random(A, rand_pool, sc_t(10));
  Kokkos::fill_random(B, rand_pool, sc_t(10));
  KokkosBlas::fill(C, one);

  // warmup
  //KokkosBlas::gemm(&ct_N, &ct_N, one, A, B, one, C);

  std::array<double, 3> perfTimes = {1e32,0.,0.}; //min, max, ave
  const auto startTime = std::chrono::high_resolution_clock::now();
  for (auto i=0; i<rep; ++i){
    Kokkos::Timer timer;

    KokkosBlas::gemm(&ct_N, &ct_N, one, A, B, one, C);

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

  const sc_t gFlopCount = double(m)*n*p*rep*2*1e-9;
  std::cout << gFlopCount/elapsed.count()  <<  " GFLOPS \n";
  double memCostMB, flopsCost = 0.;
  Complexity<sc_t>::gemm_beta_one(m, p, n, memCostMB, flopsCost);
  kokkosapp::printPerf(rep, perfTimes, memCostMB, flopsCost);
}

int main(int argc, char *argv[])
{
  // *** parse args ***
  int m = 0;
  int p = 0;
  int n = 0;
  int rep  = 0;

  if(argc>=2){
    m = std::atoi(argv[1]);
    p = std::atoi(argv[2]);
    n = std::atoi(argv[3]);
    rep  = std::atoi(argv[4]);
  }

  std::cout << m << " "
	    << p << " "
	    << n << " "
	    << rep << std::endl;

  // {
  //  gemmEigen(m, p, n, rep);
  // }

  Kokkos::initialize (argc, argv);
  gemmKokkos(m, p, n, rep);
  Kokkos::finalize();

  return 0;
}
