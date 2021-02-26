/*
//@HEADER
// ************************************************************************
//
// parser_general_section.hpp
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

#ifndef SHAXIPP_PARSER_MIXIN_GENERAL_SECTION_HPP_
#define SHAXIPP_PARSER_MIXIN_GENERAL_SECTION_HPP_

#include "yaml-cpp/parser.h"
#include "yaml-cpp/yaml.h"
#include <fstream>
#include <iostream>

template <typename scalar_t>
struct ParserGeneralSection
{
private:
  std::string meshDirName_ = "empty";
  bool checkNumDispersion_ = true;
  bool checkCfl_	   = true;
  scalar_t dt_		   = {};
  scalar_t finalTime_	   = {};
  std::size_t NSteps_		   = {};
  bool exploitForcingSparsity_ = true;

public:
  auto getMeshDir() const{ return meshDirName_; }
  auto checkDispersion() const{ return checkNumDispersion_; }
  auto checkCfl() const{ return checkCfl_; }
  auto getTimeStepSize() const{ return dt_; }
  auto getNumSteps() const{ return NSteps_; }
  auto exploitForcingSparsity() const{ return exploitForcingSparsity_; }

public:
  void parseGeneral(const std::string & inputFile)
  {
    YAML::Node node0 = YAML::LoadFile(inputFile);

    // check if section is present, otherwise throw
    const auto node = node0["general"];
    if (node)
    {
      auto entry = "meshDir";
      if (node[entry]) meshDirName_ = node[entry].as<std::string>();
      else throw std::runtime_error("Empty meshDir");

      entry = "checkNumericalDispersion";
      if (node[entry]) checkNumDispersion_ = node[entry].as<bool>();

      entry = "checkCfl";
      if (node[entry]) checkCfl_ = node[entry].as<bool>();

      entry = "dt";
      if (node[entry]) dt_ = node[entry].as<scalar_t>();
      else throw std::runtime_error("Empty dt");

      entry = "finalTime";
      if (node[entry]) finalTime_ = node[entry].as<scalar_t>();
      else throw std::runtime_error("Empty final time");

      entry = "exploitForcingSparsity";
      if (node[entry]) exploitForcingSparsity_ = node[entry].as<bool>();
    }
    else{
      throw std::runtime_error("General section in yaml input is mandatory!");
    }

    NSteps_ = static_cast<std::size_t>(finalTime_/dt_);
    this->validate();
    this->print();
  }

private:
  void validate() const{
    if (dt_<=0.){
      throw std::runtime_error("Cannot have dt <= 0");
    }

    if (finalTime_<=0.){
      throw std::runtime_error("Cannot have finalT <= 0");
    }
  }

  void print() const{
    std::cout << "\nmeshDir = "		<< meshDirName_		<< " \n"
	      << "checkDispersion = "	<< std::boolalpha << checkNumDispersion_  << " \n"
	      << "timeStep = "		<< dt_			<< " \n"
	      << "finalT = "		<< finalTime_		<< " \n"
	      << "numSteps = "		<< NSteps_		<< " \n"
	      << "exploitForcingSparsity " << exploitForcingSparsity_ << " \n";
  }
};

#endif
