/*
//@HEADER
// ************************************************************************
//
// spmv.hpp
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

#ifndef UTILS_SPMV_HPP_
#define UTILS_SPMV_HPP_

template <typename sc_t>
struct Spmv
{
  //https://books.google.it/books?id=O22cDwAAQBAJ&pg=PA97&lpg=PA97&dq=gemm+flops+supercomputing+frontiers&source=bl&ots=tiCjQSEypI&sig=ACfU3U2D2zTZHLGPVXNzIFV8SBO_jhz-1w&hl=en&sa=X&ved=2ahUKEwjNieO9nsXpAhVuwqYKHaRJDZIQ6AEwAHoECAoQAQ#v=onepage&q=gemm%20flops%20supercomputing%20frontiers&f=false

  // y = beta*y + alpha*A*x
  template <typename ordinal_t>
  static void spmv(const std::size_t nnz,
		   const std::size_t nrows,
		   double & memCostMB,
		   double & flops)
  {
    const auto scsz  = sizeof(sc_t);
    const auto ordsz = sizeof(ordinal_t);

    const double matsize = 1.0*(nnz*(scsz+ordsz) + nrows*ordsz);
    const double vec_x_r = 1.*nnz*scsz;
    const double vec_y_r = 1.*nrows*scsz;
    const double vec_y_w = 1.*nrows*scsz;

    memCostMB = (matsize + vec_x_r + vec_y_r + vec_y_w)/1024./1024.;
    flops = 2.0*nnz /*Ax*/ + nrows /*y+...*/;
  }

};
#endif
