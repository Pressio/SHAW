/*
//@HEADER
// ************************************************************************
//
// parser_forcing_section.hpp
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

#ifndef SHAXIPP_PARSER_MIXIN_FORCING_SECTION_HPP_
#define SHAXIPP_PARSER_MIXIN_FORCING_SECTION_HPP_

// forward declaration
template <typename sc_t>
class Signal;

template <typename scalar_t>
class ParserForcingSection
{

  bool multiForcing_	= false;
  bool enableRank2Mode_	= false;
  int forcingSize_ = 1;

  signalKind kind_ = {};
  std::vector<scalar_t> depths_  = {}; // km;
  std::vector<scalar_t> angles_  = {}; // degrees;
  std::vector<scalar_t> periods_ = {}; // seconds;
  std::vector<scalar_t> delays_  = {}; // seconds;

public:
  bool rank2Enabled() const {
    return enableRank2Mode_;
  }

  int getForcingSize() const {
    return forcingSize_;
  }

  bool multiForcing() const {
    return multiForcing_;
  }

  auto getSourceSignalKind() const{
    return kind_;
  }

  const std::vector<scalar_t> & viewDepths() const { return depths_; }
  const std::vector<scalar_t> & viewAngles() const { return angles_; }
  const std::vector<scalar_t> & viewPeriods() const { return periods_; }
  const std::vector<scalar_t> & viewDelays() const { return delays_; }

public:
  void parseForcing(const std::string & inputFile)
  {
    YAML::Node node = YAML::LoadFile(inputFile);

    // check if section is present, otherwise throw
    const auto node1 = node["source"];
    if (node1)
    {
      const auto node2 = node1["signal"];
      if (node2)
      {
	// check if subnode for "kind" exists
	auto n = node2["kind"];
      	if (n){
	  kind_ = stringToSignalKind(n.as<std::string>());
	}
      	else{
	  throw std::runtime_error("You must set the kind of the signal");
	}

	// check if subnode for "depth" exists
	auto n1 = node2["depth"];
       	if (n1){
	  // check if it is a scalar or an array
	  if (n1.IsScalar()){
	    depths_.push_back(n1.as<scalar_t>());
	  }
	  else{
	    depths_= n1.as<std::vector<scalar_t>>();
	  }
	}
      	else{
	  throw std::runtime_error("You must set the depth of the signal");
	}

	auto n2 = node2["angle"];
       	if (n2){
	  // check if it is a scalar or an array
	  if (n2.IsScalar()){
	    angles_.push_back(n2.as<scalar_t>());
	  }
	  else{
	    angles_ = n2.as<std::vector<scalar_t>>();
	  }
	}
      	else{
	  angles_.push_back(0.);
	  std::cout << "Angle of source not specified, default==0";
	}

	auto n3 = node2["period"];
       	if (n3){
	  // check if it is a scalar or an array
	  if (n3.IsScalar()){
	    periods_.push_back(n3.as<scalar_t>());
	  }
	  else{
	    periods_ = n3.as<std::vector<scalar_t>>();
	  }
	}
      	else{
	  throw std::runtime_error("You must set the period of the signal");
	}

	auto n4 = node2["delay"];
       	if (n4){
	  // check if it is a scalar or an array
	  if (n4.IsScalar()){
	    delays_.push_back(n4.as<scalar_t>());
	  }
	  else{
	    delays_ = n4.as<std::vector<scalar_t>>();
	  }
	}
      	else{
	  throw std::runtime_error("You must set the delay of the signal");
	}

	auto n5 = node2["forcingSize"];
       	if (n5){
	  // if forcingSize is present and >=2, then rank-2 is on
	  forcingSize_ = n5.as<int>();
	  enableRank2Mode_ = true;
	}
      }
      else{
      	throw std::runtime_error
	  ("Signal node in inputfile is missing, it is mandatory!");
      }
    }
    else{
      throw std::runtime_error
	("Cannot find source section in inputfile, must be set!");
    }

    setMultiForcingFlag();
    if (multiForcing_){
      validateMultiForcing();
      printMultiForcing();
    }
    else{
      validateSingleForcing();
      printSingleForcing();
    }
  }

private:
  // set multiForcing to true if input file demands so
  void setMultiForcingFlag()
  {
    // we have a multifFOrcing run if any of the source fields is an array
    if (depths_.size() >= 2 or
	angles_.size() >= 2 or
	periods_.size() >= 2 or
	delays_.size() >= 2)
      {
	multiForcing_ = true;
      }
  }

  void validateMultiForcing() const
  {
    if (kind_ == signalKind::unknown){
      const auto ssKstr = signalKindToString(kind_);
      const auto msg = "Invalid source kind: " + ssKstr + ", maybe you mispelled it in input file?\n";
      throw std::runtime_error(msg);
    }

    for (auto & it : depths_){
      if (it<=0.){
	throw std::runtime_error("The depth of the source must be a positive number");
      }
    }

    for (auto & it : angles_){
      if (it<0.){
	throw std::runtime_error("The angle of source must be zero or a positive value");
      }
    }

    for (auto & it : periods_){
      if (it<=0.){
	throw std::runtime_error("Period for signal must be positive. ");
      }
    }

    for (auto & it : delays_){
      if (it<0.){
	throw std::runtime_error("Delay time of the signal must be >=0");
      }
    }
  }

  void validateSingleForcing() const
  {
    if (kind_ == signalKind::unknown){
      const auto ssKstr = signalKindToString(kind_);
      const auto msg = "Invalid source kind: " + ssKstr + ", maybe you mispelled it in input file?\n";
      throw std::runtime_error(msg);
    }

    if (depths_[0]<=0.){
      throw std::runtime_error("The depth of the source must be a positive number");
    }

    if (angles_[0]<0.){
      throw std::runtime_error("The angle of source must be zero or a positive value");
    }

    if (periods_[0]<=0.){
      throw std::runtime_error("Period for signal must be positive. ");
    }

    if (delays_[0]<0.){
      throw std::runtime_error("Delay time of the signal must be >=0");
    }
  }

  void printMultiForcing() const
  {}


  void printSingleForcing() const
  {
    std::cout << std::endl;
    std::cout << "signal "	 << " "
	      << "type= "	 << signalKindToString(kind_) << " "
	      << "depth[km]= "   << depths_[0]  << " "
	      << "angle[deg]= "  << angles_[0]  << " "
	      << "period[sec]= " << periods_[0] << " "
	      << "delay[sec]= "  << delays_[0] << " \n";
  }
};

#endif
