/*
//@HEADER
// ************************************************************************
//
// parser_rom_section.hpp
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

#ifndef SHAXIPP_PARSER_MIXIN_ROM_SECTION_HPP_
#define SHAXIPP_PARSER_MIXIN_ROM_SECTION_HPP_

template <typename scalar_t>
struct ParserRomSection
{
private:
  bool romOn_			= false;

  // dummy basis is turned on when the input file
  // does not provide the basis files for vp and sp
  // which means we don't care about corrctness but
  // we just want to execute scaling/performance
  // for which we don't need real numbers
  bool dummyBasis_		= false;

  // disableCompRomJacobians
  // this is used just for scaling purpose when we are
  // not interested in computing the rom jacobians
  bool disableCompRomJacobians_ = false;

  std::size_t vpRomSize_		= {};
  std::size_t spRomSize_		= {};
  std::string vpBasisFileName_	= "empty";
  std::string spBasisFileName_	= "empty";

  bool vpBinaryBasis_ = true;
  bool spBinaryBasis_ = true;

public:
  auto enableRom() const{ return romOn_; };
  auto enableRandomDummyBasis() const{ return dummyBasis_; };
  auto disableCompRomJacobians() const{ return disableCompRomJacobians_; };

  auto readBinaryBasis(dofId dof) const{
    switch(dof){
    case dofId::vp: return vpBinaryBasis_;
    case dofId::sp: return spBinaryBasis_;
    default: return true;
    }
  };

  std::size_t getRomSize(dofId dof) const{
    switch(dof){
    case dofId::vp: return vpRomSize_;
    case dofId::sp: return spRomSize_;
    default: return 0;
    }
  };

  auto getBasisFileName(dofId dof) const{
    switch(dof){
    case dofId::vp: return vpBasisFileName_;
    case dofId::sp: return spBasisFileName_;
    default: return std::string();
    }
  };

  void parseRom(const std::string inputFile)
  {
    const YAML::Node node = YAML::LoadFile(inputFile);

    // check if section is present
    const auto romNode= node["rom"];
    if (romNode)
    {
      romOn_ = true;

      // check if we disable the comp of the rom jacobins
      const auto disableRomJac = romNode["disableCompRomJacs"];
      if (disableRomJac){
	disableCompRomJacobians_ = disableRomJac.as<bool>();
      }

      // get inputs for velocity
      const auto veloNode = romNode["velocity"];
      if (veloNode){
	auto entry = "numModes";
	if (veloNode[entry]) vpRomSize_ = veloNode[entry].as<std::size_t>();
	else throw std::runtime_error("You must set # of modes for velocity");

	entry = "modesFile";
	if (veloNode[entry]) vpBasisFileName_ = veloNode[entry].as<std::string>();

	if (veloNode["binary"]) vpBinaryBasis_ = veloNode["binary"].as<bool>();
      }
      else
	throw std::runtime_error("Cannot find rom inputs for velocity");

      // get inputs for stress
      const auto stressNode = romNode["stress"];
      if (stressNode){
	auto entry = "numModes";
	if (stressNode[entry]) spRomSize_ = stressNode[entry].as<std::size_t>();
	else throw std::runtime_error("You must set # of modes for stress");

	entry = "modesFile";
	if (stressNode[entry]) spBasisFileName_ = stressNode[entry].as<std::string>();

	if (stressNode["binary"]) spBinaryBasis_ = stressNode["binary"].as<bool>();
      }
      else
	throw std::runtime_error("Cannot find rom inputs for stress");

      // if {vp,sp}BasisFileName are empty, it means we use dummy basis
      if (vpBasisFileName_ == "empty" or spBasisFileName_ == "empty"){
	std::cout << "You did not set the basisfiles so I am using dummy basis" << std::endl;
	dummyBasis_ = true;
	}

      this->validate();
    }

    this->print();
  }

private:
  void validate() const
  {
    if (vpRomSize_<=0) throw std::runtime_error("Cannot have vpRomSize <=0 ");
    if (spRomSize_<=0) throw std::runtime_error("Cannot have spRomSize <=0 ");
  }

  void print() const
  {
    if (romOn_){
      std::cout << "vpRomSize_ = "		<< vpRomSize_		<< " \n"
		<< "spRomSize_ = "		<< spRomSize_		<< " \n"
		<< "vpBasisFileName_ = "	<< vpBasisFileName_	<< " \n"
		<< "spBasisFileName_ = "	<< spBasisFileName_	<< "\n"
		<< "---------------------------------------------------\n";
    }
  }
};

#endif
