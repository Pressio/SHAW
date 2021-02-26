/*
//@HEADER
// ************************************************************************
//
// material_model_create.hpp
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

#ifndef MATERIAL_MODEL_CREATE_HPP_
#define MATERIAL_MODEL_CREATE_HPP_

template<typename scalar_t, typename parser_t, typename mesh_info_t>
std::shared_ptr<MaterialModelBase<scalar_t>>
createMaterialModel(const parser_t & parser, const mesh_info_t & meshInfo)
{
  auto kind = parser.getMaterialModelKind();
  if (kind == materialModelKind::unilayer){
    using ret_t = UnilayerMaterialModel<scalar_t, parser_t>;
    return std::make_shared<ret_t>(parser, meshInfo);
  }
  else if (kind == materialModelKind::bilayer){
    using ret_t = BilayerMaterialModel<scalar_t, parser_t>;
    return std::make_shared<ret_t>(parser, meshInfo);
  }
  else if (kind == materialModelKind::prem){
    using ret_t = PremMaterialModel<scalar_t, parser_t>;
    return std::make_shared<ret_t>(parser, meshInfo);
  }
  else
    throw std::runtime_error("Cannot create material object, invalid material model kind");
}

#endif
