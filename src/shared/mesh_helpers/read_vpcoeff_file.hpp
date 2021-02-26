/*
//@HEADER
// ************************************************************************
//
// read_vpcoeff_file.hpp
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

#ifndef READ_VPCOEFF_FILE_HPP_
#define READ_VPCOEFF_FILE_HPP_

template <typename sc_t, typename coeffs_t>
void readFullMeshCoeffFile(std::string meshDir, dofId dofid, coeffs_t & coeffs)
{
  const std::string dofName = dofIdToString(dofid);
  const std::string filePath = meshDir + "/coeff_" + dofName + ".dat";

  std::cout << "Reading Vp stencil coeffs for " << dofName << "...";

  std::ifstream foundFile(filePath);
  if(!foundFile){
    throw std::runtime_error(" file: " + filePath + "not found");
  }

  std::ifstream source;
  source.open( filePath, std::ios_base::in);
  std::string line, colVal;
  while (std::getline(source, line) ){
    std::istringstream ss(line);
    // first col contains the gid
    ss >> colVal; const auto currVpGid = std::stoi(colVal);
    // store the coeffcients
    for (auto i=0; i<coeffs.extent(1); ++i){
      ss >> colVal;
      coeffs(currVpGid, i) = std::stod(colVal);
    }
  }//while

  source.close();
  std::cout << "Done" << std::endl;
}//end readCoeffsVp

#endif
