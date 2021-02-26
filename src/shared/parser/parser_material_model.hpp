/*
//@HEADER
// ************************************************************************
//
// parser_material_model.hpp
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

#ifndef SHAXIPP_PARSER_MIXIN_MATERIAL_MODEL_HPP_
#define SHAXIPP_PARSER_MIXIN_MATERIAL_MODEL_HPP_

template <typename scalar_t>
struct ParserMaterialModel
{
  /* we use:
   * 2. 1d array to store the nominal locations of the discontinuities
   * 1. 2d array to store the parametrization of the material properties
   *
   * The 2d array is such that each row identifies a layer,
   * and each col contains the paramtrization coefficients for a polyn.
   * For the time being, we limit the parametrization to be at most 2nd-order.
   * THis means we can have at most 3 coefficients.
   *
   * The 1d array with discotinuities is such that each row contains
   * the depth (km) of the discontinuity.
   *
   * For example, let's say that we have a single layer model, which yields:
   * density_array = [c0, c1, c2] such that rho(depth) = c0 + depth*c1 + c2*depth^2
   * vs_array	   = [d0, d1, d2] such that vs(depth)  = d0 + depth*d1 + d2*depth^2
   * discon_array  = [0] because the discontinuity here is just the surface
   *
   * For example, let's say that we have a two-layer model, with the
   * discontinuity located at 660 km of depth, this yields:
   * density_array = [c0, c1, c2] such that rho(depth) = c0 + depth*c1 + c2*depth^2
   * vs_array	   = [d0, d1, d2] such that vs(depth)  = d0 + depth*d1 + d2*depth^2
   * discon_array  = [660]
   *
   */

  using poly_coeff_t     = std::array<scalar_t, 3>;
  using profile_params_t = std::vector<poly_coeff_t>;
  using discont_depth_t  = std::vector<scalar_t>;

private:
  materialModelKind matModKind_	= materialModelKind::unknown;

  profile_params_t density_ = {};
  profile_params_t vs_	    = {};
  discont_depth_t discont_  = {};

public:
  auto getMaterialModelKind() const{ return matModKind_; }

  const discont_depth_t  & viewDiscontinuityDepthsKm() const{ return discont_; }
  const profile_params_t & viewDensityParametrization() const{ return density_; }
  const profile_params_t & viewVelocityParametrization() const{ return vs_; }

  void parseMaterial(const std::string & inputFile)
  {
    const YAML::Node node0 = YAML::LoadFile(inputFile);

    // check if present, otherwise throw
    if (node0["material"])
    {
      const auto node = node0["material"];

      const auto matKindStr = node["kind"].as<std::string>();
      if (matKindStr == "unilayer" or
      	  matKindStr == "Unilayer" or
      	  matKindStr == "UniLayer")
      {
      	matModKind_ = materialModelKind::unilayer;
	density_.resize(1); vs_.resize(1); discont_.resize(1);
	this->parseUnilayer(node);
	this->validateSingleLayer();
      }

      else if (matKindStr == "bilayer" or
      	       matKindStr == "Bilayer" or
      	       matKindStr == "BiLayer")
      {
      	matModKind_ = materialModelKind::bilayer;
	density_.resize(2); vs_.resize(2); discont_.resize(2);
	this->parseBilayer(node);
	this->validateBiLayer();
      }

      else if (matKindStr == "prem" or matKindStr == "PREM"){
      	matModKind_ = materialModelKind::prem;
      }

      else if (matKindStr == "custom"){
      	matModKind_ = materialModelKind::custom;
      }

      else{
      	throw std::runtime_error("Unknown material model kind, maybe mispelled it?");
      }
    }
    else{
      throw std::runtime_error("Cannot find material section in yaml");
    }

    this->print();
  }

private:
  void parseUnilayer(const YAML::Node & node)
  {
    // for a single layer, there is no depth entry since it is a layer below surface
    discont_[0] = 0.;

    // get density coeff
    auto densC = node["layer"]["density"].as<std::vector<scalar_t>>();
    for (std::size_t i=0; i<densC.size(); ++i) density_[0][i] = densC[i];
    // get vs coeff
    auto vsC = node["layer"]["velocity"].as<std::vector<scalar_t>>();
    for (std::size_t i=0; i<vsC.size(); ++i) vs_[0][i] = vsC[i];
  }

  void parseBilayer(const YAML::Node & node)
  {
    // get data for first layer
    // (this does not have depth entry since it is layer below surface)
    discont_[0] = 0.;

    auto densC = node["layer1"]["density"].as<std::vector<scalar_t>>();
    for (std::size_t i=0; i<densC.size(); ++i) density_[0][i] = densC[i];

    auto vsC = node["layer1"]["velocity"].as<std::vector<scalar_t>>();
    for (std::size_t i=0; i<vsC.size(); ++i) vs_[0][i] = vsC[i];

    // get data for second layer (needs to have a depth entry to identify discont)
    discont_[1] = node["layer2"]["depth"].as<scalar_t>();

    densC = node["layer2"]["density"].as<std::vector<scalar_t>>();
    for (std::size_t i=0; i<densC.size(); ++i) density_[1][i] = densC[i];

    vsC = node["layer2"]["velocity"].as<std::vector<scalar_t>>();
    for (std::size_t i=0; i<vsC.size(); ++i) vs_[1][i] = vsC[i];
  }

  void print() const
  {
    std::cout << "\n";
    std::cout << "matModelType = "	<< materialModelKindToString(matModKind_) << " \n";
    if (matModKind_==materialModelKind::unilayer){
      std::cout << "rho(c0,c1,c2) = "
		<< density_[0][0] << " "
		<< density_[0][1] << " "
		<< density_[0][2] << "\n";
      std::cout << "vs(c0,c1,c2)  = "
		<< vs_[0][0] << " "
		<< vs_[0][1] << " "
		<< vs_[0][2] << "\n";
    }

    if (matModKind_==materialModelKind::bilayer){
      std::cout << std::endl;
      std::cout << "Layer 1 \n";
      std::cout << "rho(c0,c1,c2) = "
		<< density_[0][0] << " "
		<< density_[0][1] << " "
		<< density_[0][2] << "\n";
      std::cout << "vs(c0,c1,c2)  = "
		<< vs_[0][0] << " "
		<< vs_[0][1] << " "
		<< vs_[0][2] << "\n";
      std::cout << "Layer 2 \n";
      std::cout << "Starting depth [km] " << discont_[1] << std::endl;
      std::cout << "rho(c0,c1,c2) = "
		<< density_[1][0] << " "
		<< density_[1][1] << " "
		<< density_[1][2] << "\n";
      std::cout << "vs(c0,c1,c2)  = "
		<< vs_[1][0] << " "
		<< vs_[1][1] << " "
		<< vs_[1][2] << "\n";
    }

    if (matModKind_==materialModelKind::prem) std::cout << "Using PREM model" << std::endl;

    if (matModKind_==materialModelKind::custom) std::cout << "Using custom model" << std::endl;
  }

  void validateSingleLayer() const
  {
    //tbd
  }

  void validateBiLayer() const
  {
    //TBD
  }

};

#endif
