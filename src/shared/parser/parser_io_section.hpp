/*
//@HEADER
// ************************************************************************
//
// parser_io_section.hpp
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

#ifndef SHAXIPP_PARSER_MIXIN_IO_SECTION_HPP_
#define SHAXIPP_PARSER_MIXIN_IO_SECTION_HPP_

template <typename scalar_t>
struct ParserIoSection
{
  using receivers_loc_t = std::vector<scalar_t>;

private:
  enum class writeMode{ binary, ascii };

  std::string writeModeToString(writeMode e) const{
    switch (e){
    default:	     return "unknown";
    case writeMode::binary: return "binary";
    case writeMode::ascii:  return "ascii";
    }
  }

  // *** snapshot matrix ***
  bool enableSnapMatrix_      = false;
  writeMode snapWriteMode_    = writeMode::binary;
  std::size_t vpSnapFreq_     = 0;
  std::size_t spSnapFreq_     = 0;
  std::string vpSnapFileName_ = "snaps_vp";
  std::string spSnapFileName_ = "snaps_sp";

  // *** seismogram ***
  bool enableSeismo_		  = false;
  writeMode seismoWriteMode_      = writeMode::binary;
  std::size_t seismoFreq_	  = 0;
  std::string seismogramFileName_ = "seismogram";
  receivers_loc_t receiversLocs_  = {5,30,60,90,120,150,175};

public:
  auto enableSnapshotMatrix() const{ return enableSnapMatrix_; }
  auto writeSnapshotsBinary()  const{ return snapWriteMode_ == writeMode::binary; }

  std::size_t getSnapshotFreq(const dofId & dof) const{
    switch(dof){
    case dofId::vp: return vpSnapFreq_;
    case dofId::sp: return spSnapFreq_;
    default: return 0;
    }
  }

  auto getSnapshotFileName(const dofId & dof) const{
    switch(dof){
    case dofId::vp: return vpSnapFileName_;
    case dofId::sp: return spSnapFileName_;
    default: return std::string();
    }
  }

  auto enableSeismogram()     const{ return enableSeismo_; }
  auto writeSeismogramBinary() const{ return seismoWriteMode_ == writeMode::binary; }
  auto getSeismogramFileName() const{ return seismogramFileName_; }
  auto getSeismoFreq() const{ return seismoFreq_; }
  auto getSeismoReceiversAnglesDeg() const{ return receiversLocs_; }

public:
  void parseIo(const std::string & inputFile)
  {
    const YAML::Node node = YAML::LoadFile(inputFile);

    // check if io section is present
    // io is optional so no need to catch the else
    const auto ioNode= node["io"];
    if (ioNode)
    {
      // check if snapshotMatrix entry node is present
      const auto snapMatNode = ioNode["snapshotMatrix"];
      if (snapMatNode){
	enableSnapMatrix_ = true;
	this->parseSnapshotMatrixInputs(snapMatNode);
      }

      // seismo entry node
      const auto seismoNode = ioNode["seismogram"];
      if (seismoNode){
	enableSeismo_ = true;
	this->parseSeismoInputs(seismoNode);
      }

      this->validate();
    }

    this->print();
  }

private:
  void parseSnapshotMatrixInputs(const YAML::Node & node)
  {
    if (node["binary"]){
      auto useBinary = node["binary"].as<bool>();
      snapWriteMode_ = useBinary ? writeMode::binary : writeMode::ascii;
    }

    const auto veloNode = node["velocity"];
    if (veloNode)
    {
      auto entry = "freq";
      if (veloNode[entry]){
	vpSnapFreq_ = veloNode[entry].as<std::size_t>();
      }
      else{
	throw std::runtime_error("You must set snapshot freq for velocity");
      }

      entry = "fileName";
      if (veloNode[entry]) vpSnapFileName_ = veloNode[entry].as<std::string>();
    }
    else{
      throw std::runtime_error("Cannot find io inputs for velocity");
    }

    const auto stressNode = node["stress"];
    if (stressNode)
    {
      auto entry = "freq";
      if (stressNode[entry]){
	spSnapFreq_ = stressNode[entry].as<std::size_t>();
      }
      else{
	throw std::runtime_error("You must set snapshot freq for stress");
      }

      entry = "fileName";
      if (stressNode[entry]) spSnapFileName_ = stressNode[entry].as<std::string>();
    }
    else
      throw std::runtime_error("Cannot find io inputs for stress");
  }

  void parseSeismoInputs(const YAML::Node & node)
  {
    if (node["binary"]){
      auto useBinary = node["binary"].as<bool>();
      seismoWriteMode_ = useBinary ? writeMode::binary : writeMode::ascii;
    }

    if (node["fileName"]){
      seismogramFileName_ = node["fileName"].as<std::string>();
    }

    auto entry = "freq";
    if (node[entry]) {
      seismoFreq_ = node[entry].as<std::size_t>();
    }
    else{
      throw std::runtime_error("You must set freq for seismogram");
    }

    // check for receivers
    const auto recNode = node["receivers"];
    if (recNode){
      // check if we have a vector or a file
      if (recNode.IsSequence()){
	receiversLocs_ = recNode.as<receivers_loc_t>();
      }
    }
  }

  void validate() const
  {
    if (enableSnapMatrix_){
      if (vpSnapFreq_<=0) throw std::runtime_error("cannot have vpSnapshotsFreq <=0 ");
      if (spSnapFreq_<=0) throw std::runtime_error("cannot have spSnapshotsFreq <=0 ");
    }

    if (enableSeismo_){
      if (seismoFreq_<=0) throw std::runtime_error("cannot have seismoFreq <=0 ");
    }
  }

  void print() const
  {
    std::cout << std::endl;
    std::cout << "enableSnapshotMatrix = " << std::boolalpha << enableSnapMatrix_ << " \n";
    if (enableSnapMatrix_){
      std::cout << "mode = " << writeModeToString(snapWriteMode_)  << " \n"
		<< "vpSnapshotsFileName_ = "  << vpSnapFileName_  << " \n"
		<< "spSnapshotsFileName_ = "  << spSnapFileName_  << " \n"
		<< "vpSnapshotsFreq_ = "      << vpSnapFreq_	  << " \n"
		<< "spSnapshotsFreq_ = "      << spSnapFreq_	  << " \n";
    }

    std::cout << "enableSeimogram = " << std::boolalpha << enableSeismo_ << " \n";
    if (enableSeismo_){
      std::cout << "mode = "	<< writeModeToString(seismoWriteMode_) << " \n";
      std::cout << "Freq_ = "	<< seismoFreq_			      << " \n";
      std::cout << "Locations: ";
      for (const auto & it : receiversLocs_) std::cout << it << " ";
    }
  }

};
#endif
