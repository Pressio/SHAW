
#ifndef SHAXIPP_KOKKOS_PRINT_PERF_HPP_
#define SHAXIPP_KOKKOS_PRINT_PERF_HPP_

#include <array>
#include <iostream>

namespace kokkosapp{

void printPerf(std::size_t loops,
	       const std::array<double,3> perfTimes, // times are in seconds
	       const double memCostMB,
	       const double flops)
{
  const auto minTimeSec	  = perfTimes[0];
  const auto maxTimeSec   = perfTimes[1];
  const auto totalTimeSec = perfTimes[2];
  const auto aveTimeSec   = totalTimeSec/(double)loops;
  const double memCostGB  = memCostMB/1024.;

  std::cout << "flops = " << flops << std::endl;
  std::cout << "memMB = " << memCostMB << std::endl;
  std::cout << "flops/bytes = " << flops/(memCostMB*1024.*1024.) << std::endl;

  printf("aveBandwidth(GB/s) = %8.2lf \n", memCostGB/aveTimeSec);
  printf("minBandwidth(GB/s) = %8.2lf \n", memCostGB/maxTimeSec);
  printf("maxBandwidth(GB/s) = %8.2lf \n", memCostGB/minTimeSec);

  printf("aveGFlop = %7.4lf \n", flops/aveTimeSec/1e9);
  printf("minGFlop = %7.4lf \n", flops/maxTimeSec/1e9);
  printf("maxGFlop = %7.4lf \n", flops/minTimeSec/1e9);

  printf("totTime(se) = %7.4lf \n", perfTimes[2]);
  printf("aveTime(ms) = %7.4lf \n", aveTimeSec*1000);
  printf("maxTime(ms) = %7.4lf \n", maxTimeSec*1000);
  printf("minTime(ms) = %7.4lf \n",  minTimeSec*1000);
}

}// end namespace
#endif
