/*
//@HEADER
// ************************************************************************
//
// axpy.hpp
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

#ifndef UTILS_AXPY_HPP_
#define UTILS_AXPY_HPP_

template <typename sc_t>
struct Axpy
{
  // y = y + alpha*x

  // y,x rank-2 R^(m,n) and alpha = rank-1 R^n
  static void axpy(const std::size_t m,
		   const std::size_t n,
		   double & memCostMB,
		   double & flops)
  {
    const auto scsz  = sizeof(sc_t);
    const double reads = 1.0*(2.*m*n + 1.*n)*scsz;
    const double write = 1.0*(m*n)*scsz;
    memCostMB = (reads+write)/1024./1024.;
    flops = 2.*n;
  }

  // alpha = scalar
  static void axpy(const std::size_t n,
		   double & memCostMB,
		   double & flops)
  {
    const auto scsz  = sizeof(sc_t);
    memCostMB = (3.0*n*scsz)/1024./1024.;
    flops = 2.*n;
  }

  // alpha = 0
  static void axpy_alpha_zero(const std::size_t n,
			      double & memCostMB,
			      double & flops)
  {
    memCostMB = 0.;
    flops = 0.;
  }

  // alpha=1
  static void axpy_alpha_one(const std::size_t n,
			     double & memCostMB,
			     double & flops)
  {
    const auto scsz  = sizeof(sc_t);
    memCostMB = (3.0*n*scsz)/1024./1024.;
    flops = 1.*n;
  }

};
#endif
