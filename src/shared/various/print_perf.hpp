/*
//@HEADER
// ************************************************************************
//
// print_perf.hpp
//                     		Pressio/SHAW
//                         Copyright 2019
// National Technology & Engineering Solutions of Sandia, LLC (NTESS)
//
// Under the terms of Contract DE-NA0003525 with NTESS, the
// U.S. Government retains certain rights in this software.
//
// Pressio is licensed under BSD-3-Clause terms of use:
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
// contributors may be used to endorse or promote products derived
// from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Francesco Rizzi (fnrizzi@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

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
