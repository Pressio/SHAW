/*
//@HEADER
// ************************************************************************
//
// read_graph_file.hpp
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

#ifndef READ_GRAPH_FILE_HPP_
#define READ_GRAPH_FILE_HPP_

namespace{

template <typename sc_t, typename graph_t, typename coords_t, typename cot_t, typename labels_t>
void _readFullMeshGraphFileImpl(const dofId dofid,
			       const std::string meshDir,
			       graph_t & graph,
			       coords_t & coords,
			       cot_t & cot,
			       labels_t & labels,
			       bool readLabels)
{
  constexpr auto zero	  = constants<sc_t>::zero();
  constexpr auto one	  = constants<sc_t>::one();
  constexpr auto thousand= constants<sc_t>::thousand();

  const std::string dofName = dofIdToString(dofid);
  std::cout << "Reading graph for " << dofName << " ...";

  const std::string filePath = meshDir + "/graph_" + dofName + ".dat";
  std::ifstream foundFile(filePath);
  if(!foundFile){
    const std::string errMsg = "graph file: " + filePath + " not found";
    throw std::runtime_error(errMsg);
  }

  std::ifstream source;
  source.open( filePath, std::ios_base::in);
  std::string line, colVal;
  while (std::getline(source, line) ){
    std::istringstream ss(line);

    // first col contains the gid of current point
    ss >> colVal; const auto currGid = std::stoi(colVal);
    // store it as first entry of the adjecency list
    graph(currGid, 0) = currGid;

    if (readLabels){
      ss >> colVal;
      labels(currGid) = std::stoi(colVal);
    }

    // read col indicating if a vp or srp point is on symmetry axes
    ss >> colVal;
    const int onSymAxis = std::stoi(colVal);

    // *** store the theta coordinate of this point ***
    ss >> colVal;
    coords(currGid,0) = std::stod(colVal);

    // *** store the radius of this point ***
    // multiply by 1000 to convert radius from km to m
    ss >> colVal;
    coords(currGid, 1) = one/(std::stod(colVal)*thousand);

    // *** compute cotangent for current point ***
    if (onSymAxis==1)
      cot(currGid) = zero;
    else
      cot(currGid) = computeCotangent(coords(currGid,0));

    // When we deal with Vp, all connected points contain stress dofs.
    // When we deal with Sp, all connected points contain velocity dofs.
    for (auto i=1; i<=graph.extent(1)-1; ++i){
      // get the gid for the connected node
      ss >> colVal; const auto testGid = std::stoi(colVal);
      // add gid of connected node
      graph(currGid, i) = testGid;
    }
  }//while

  source.close();
  std::cout << "Done" << std::endl;
}
}//anonym namespace

template <typename sc_t, typename graph_t, typename coords_t, typename cot_t>
void readFullMeshGraphFile(const std::string meshDir, dofId dofid,
			   graph_t & graph, coords_t & coords, cot_t & cot)
{
  using labels_t = cot_t;
  labels_t dummy;
  _readFullMeshGraphFileImpl<sc_t>(dofid, meshDir, graph, coords, cot, dummy, false);
}

template <typename sc_t, typename graph_t, typename coords_t, typename cot_t, typename labels_t>
void readFullMeshGraphFile(const std::string meshDir, dofId dofid,
			   graph_t & graph, coords_t & coords, cot_t & cot, labels_t & labels)
{
  _readFullMeshGraphFileImpl<sc_t>(dofid, meshDir, graph, coords, cot, labels, true);
}

#endif
