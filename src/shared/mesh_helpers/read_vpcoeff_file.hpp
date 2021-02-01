
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
