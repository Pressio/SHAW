
#ifndef READ_GRAPH_FILE_HPP_
#define READ_GRAPH_FILE_HPP_

namespace{

// template <typename sc_t, typename graph_t, typename vector_t, typename labels_t>
// void readGraphFileImpl(const std::string dofName, // vp or sp
// 		       const std::string meshDir,
// 		       graph_t & graph,
// 		       vector_t & thVp, vector_t & rrVpInv, vector_t & cotVp,
// 		       vector_t & thSp, vector_t & rrSpInv, vector_t & cotSp,
// 		       labels_t & labels)
// {
//   constexpr auto zero	  = constants<sc_t>::zero();
//   constexpr auto one	  = constants<sc_t>::one();
//   constexpr auto thousand= constants<sc_t>::thousand();

//   if (dofName != "vp" and dofName != "sp")
//     throw std::runtime_error("Invalid dofName passed to reagGraphFileImpl");

//   std::cout << "Reading graph for " << dofName << std::endl;

//   const std::string filePath = meshDir + "/graph_" + dofName + ".dat";
//   std::ifstream foundFile(filePath);
//   if(!foundFile){
//     std::cout << "graph file: " << filePath << " not found" << std::endl;
//     exit(EXIT_FAILURE);
//   }

//   std::ifstream source;
//   source.open( filePath, std::ios_base::in);
//   std::string line;
//   int64_t count = -1;
//   std::string colVal;
//   while (std::getline(source, line) ){
//     std::istringstream ss(line);

//     // first col contains the gid of current point
//     ss >> colVal; const auto currGid = std::stoi(colVal);
//     // store it as first entry of the adjecency list
//     graph(++count, 0) = currGid;

//     // second column of file has the label, 0=vp, 1=srp, 2=stp
//     ss >> colVal; labels(count) = std::stoi(colVal);

//     // *** store the theta coordinate of this point ***
//     ss >> colVal;
//     if (dofName == "vp")
//       thVp(currGid) = std::stod(colVal);
//     else
//       thSp(currGid) = std::stod(colVal);

//     // *** store the radius of this point ***
//     // multiply by 1000 to convert radius from km to m
//     ss >> colVal;
//     if (dofName == "vp")
//       rrVpInv(currGid) = one/(std::stod(colVal)*thousand);
//     else
//       rrSpInv(currGid) = one/(std::stod(colVal)*thousand);

//     // *** compute cotangent for current point ***
//     if (dofName == "vp"){
//       cotVp(currGid) = computeCotangent(thVp(currGid));
//     }
//     else{
//       cotSp(currGid) = computeCotangent(thSp(currGid));
//     }

//     // store the adjecency list of the current point
//     int numNeigh = 4;
//     if (dofName=="sp")
//       numNeigh = 2;

//     // When we deal with Vp, all connected points contain stress dofs.
//     // When we deal with Sp, all connected points contain velocity dofs.
//     for (auto i=1; i<=numNeigh; ++i){
//       // get the gid for the connected node
//       ss >> colVal; const auto testGid = std::stoi(colVal);
//       // add gid of connected node
//       graph(count, i) = testGid;

//       // get theta for the connected node
//       ss >> colVal;
//       if (dofName == "vp")
//       	thSp(testGid) = std::stod(colVal);
//       else
//       	thVp(testGid) = std::stod(colVal);

//       //compute cotangent
//       if (dofName == "vp"){
//       	cotSp(testGid) = computeCotangent(thSp(testGid));
//       }
//       else{
//       	cotVp(testGid) = computeCotangent(thVp(testGid));
//       }

//       // get radius for the connected node
//       ss >> colVal;
//       if (dofName == "vp")
//       	rrSpInv(testGid) = one/(std::stod(colVal)*thousand);
//       else
//       	rrVpInv(testGid) = one/(std::stod(colVal)*thousand);
//     }
//   }//while

//   source.close();
//   std::cout << "End Reading graph for " << dofName << std::endl;
// }//end readGraphVp


template <typename sc_t, typename graph_t, typename coords_t, typename cot_t, typename labels_t>
void readFullMeshGraphFileImpl(const dofId dofid,
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
}//end readGraph

}//anonym namespace


template <typename sc_t, typename graph_t, typename coords_t, typename cot_t>
void readFullMeshGraphFile(const std::string meshDir, dofId dofid,
			   graph_t & graph, coords_t & coords, cot_t & cot)
{
  using labels_t = cot_t;
  labels_t dummy;
  readFullMeshGraphFileImpl<sc_t>(dofid, meshDir, graph, coords, cot, dummy, false);
}

template <typename sc_t, typename graph_t, typename coords_t, typename cot_t, typename labels_t>
void readFullMeshGraphFile(const std::string meshDir, dofId dofid,
			   graph_t & graph, coords_t & coords, cot_t & cot, labels_t & labels)
{
  readFullMeshGraphFileImpl<sc_t>(dofid, meshDir, graph, coords, cot, labels, true);
}


// template <typename sc_t, typename graph_t, typename vector_t, typename labels_t>
// void readVpGraphFile(const std::string meshDir,
// 		     graph_t & graph,
// 		     vector_t & thVp, vector_t & rrVpInv, vector_t & cotVp,
// 		     vector_t & thSp, vector_t & rrSpInv, vector_t & cotSp,
// 		     labels_t & labels)
// {
//   readGraphFileImpl<sc_t>("vp", meshDir, graph, thVp, rrVpInv, cotVp, thSp, rrSpInv, cotSp, labels);
// }


// template <typename sc_t, typename graph_t, typename vector_t, typename labels_t>
// void readSpGraphFile(const std::string meshDir,
// 		     graph_t & graph,
// 		     vector_t & thVp, vector_t & rrVpInv, vector_t & cotVp,
// 		     vector_t & thSp, vector_t & rrSpInv, vector_t & cotSp,
// 		     labels_t & labels)
// {
//   readGraphFileImpl<sc_t>("sp", meshDir, graph, thVp, rrVpInv, cotVp, thSp, rrSpInv, cotSp, labels);
// }

#endif
