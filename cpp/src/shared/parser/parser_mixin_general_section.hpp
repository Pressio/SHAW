
#ifndef SHAXIPP_PARSER_MIXIN_GENERAL_SECTION_HPP_
#define SHAXIPP_PARSER_MIXIN_GENERAL_SECTION_HPP_

#include "yaml-cpp/parser.h"
#include "yaml-cpp/yaml.h"
#include <fstream>
#include <iostream>

template <typename scalar_t, typename int_t>
struct ParserGeneralSectionMixin
{
private:
  std::string meshDirName_ = "empty";
  bool checkNumDispersion_ = true;
  bool checkCfl_	   = true;
  scalar_t dt_		   = {};
  scalar_t finalTime_	   = {};
  int_t NSteps_		   = {};
  bool includeMatPropInJacobian_ = true;
  bool exploitForcingSparsity_ = true;

public:
  auto getMeshDir() const{ return meshDirName_; }
  auto checkDispersion() const{ return checkNumDispersion_; }
  auto checkCfl() const{ return checkCfl_; }
  auto getTimeStepSize() const{ return dt_; }
  auto getNumSteps() const{ return NSteps_; }
  auto includeMatPropInJac() const{ return includeMatPropInJacobian_; }
  auto exploitForcingSparsity() const{ return exploitForcingSparsity_; }

public:
  void parseGeneral(const std::string inputFile)
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

      entry = "includeMatPropInJacobian";
      if (node[entry]) includeMatPropInJacobian_ = node[entry].as<bool>();

      entry = "exploitForcingSparsity";
      if (node[entry]) exploitForcingSparsity_ = node[entry].as<bool>();

      if (!includeMatPropInJacobian_ and exploitForcingSparsity_){
	std::cout << "WARNING: When includeMatPropInJacobian_ == false, \
I cannot exploit in a suitable way the sparsity of forcing, \
so setting exploitForcingSparsity_ = false \n";
      }

    }
    else
      throw std::runtime_error("Cannot find general section in yaml");

    NSteps_ = static_cast<int_t>(finalTime_/dt_);
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
	      << "includeMatPropInJacobian " << includeMatPropInJacobian_ << " \n"
	      << "exploitForcingSparsity " << exploitForcingSparsity_ << " \n";
  }
};

#endif
