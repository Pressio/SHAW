/*
//@HEADER
// ************************************************************************
//
// gemv.hpp
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

#ifndef UTILS_GEMV_HPP_
#define UTILS_GEMV_HPP_

template <typename sc_t>
struct Gemv
{
  // y = beta*y + alpha*A*x
  // y in R^m
  // A in R^(m,n)
  // x in R^n

  // general case: y = beta*y + alpha*op(A)*x
  static void gemv(const std::size_t m,
		   const std::size_t n,
		   double & memCostMB,
		   double & flops)
  {
    const auto scsz  = sizeof(sc_t);

    const double matread = 1.0*m*n*scsz;
    const double vec_r   = 1.*(m+n)*scsz;
    const double vec_w   = 1.*m*scsz;
    memCostMB = (matread + vec_r + vec_w)/1024./1024.;
    flops = 2.*m*n + 2.*m;
  }

  // y = y + alpha*op(A)*x
  static void gemv_beta_one(const std::size_t m,
			    const std::size_t n,
			    double & memCostMB,
			    double & flops)
  {
    const auto scsz  = sizeof(sc_t);
    const double matread = 1.0*m*n*scsz;
    const double vec_r   = 1.*(m+n)*scsz;
    const double vec_w   = 1.*m*scsz;
    memCostMB = (matread + vec_r + vec_w)/1024./1024.;
    flops = 2.*m*n + 1.*m;
  }

};
#endif
