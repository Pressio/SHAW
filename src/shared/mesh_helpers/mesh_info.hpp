/*
//@HEADER
// ************************************************************************
//
// mesh_info.hpp
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

#ifndef READ_MESH_FILE_INFO_HPP_
#define READ_MESH_FILE_INFO_HPP_

template <typename sc_t, bool isFullMesh = true>
class MeshInfo;

// specialize for when we deal with a full mesh
template <typename sc_t>
class MeshInfo<sc_t, true>
{
public:
  using ordinal_type = std::size_t;
  using domain_bounds_type = std::array<sc_t, 4>;

private:
  // dir where to find mesh files
  std::string meshDir_ = {};

  // domain bounds: theta_left (deg), theta_right (deg), r_cmb (m), r_surf (m)
  domain_bounds_type domainBounds_ = {};

  // spacing in theta (rad) and r (m) direction
  sc_t dth_{};
  sc_t drr_{};

  // inverse spacing in theta (rad) and r (m) direction
  sc_t dthInv_{};
  sc_t drrInv_{};

  // numGptVp = number of points for the velocity (Vp)
  ordinal_type numGptVp_ = {};

  // number of stress points
  ordinal_type numGptSp_ = {};

  // num of points along R/theta
  ordinal_type numPtsAlongR_ = {};
  ordinal_type numPtsAlongTh_ = {};

public:
  MeshInfo(const std::string & meshDir)
  : meshDir_{meshDir}{
    this->readMeshInfoFile();
  }

  const std::string  getMeshDir() const{ return meshDir_; }

  const sc_t  getAngularSpacing() const{ return dth_; }
  const sc_t  getRadialSpacing() const{ return drr_; }
  const sc_t  getAngularSpacingInverse() const{ return dthInv_; }
  const sc_t  getRadialSpacingInverse() const{ return drrInv_; }

  const sc_t  getMinRadius() const{ return domainBounds_[2];}
  const sc_t  getMinRadiusKm() const{ return domainBounds_[2]/1000.;}
  const sc_t  getMaxRadius() const{ return domainBounds_[3];}
  const sc_t  getMaxRadiusKm() const{ return domainBounds_[3]/1000.;}

  const sc_t  getMaxArc() const{ return dth_ * domainBounds_[3]; }
  const sc_t  getMinArc() const{ return dth_ * domainBounds_[2]; }
  const domain_bounds_type & viewDomainBounds() const{ return domainBounds_; }

  ordinal_type getNumPtsAlongR() const{ return numPtsAlongR_; }
  ordinal_type getNumPtsAlongTheta() const{ return numPtsAlongTh_; }
  ordinal_type getNumVpPts() const{ return numGptVp_; }
  ordinal_type getNumSpPts() const{ return numGptSp_; }

private:
  void readMeshInfoFile()
  {
    constexpr auto one   = constants<sc_t>::one();
    constexpr auto thousand= constants<sc_t>::thousand();

    const std::string filePath = meshDir_ + "/mesh_info.dat";
    std::ifstream foundFile(filePath);
    if(!foundFile){
      std::cout << "File named: mesh_info.dat not found" << std::endl;
      exit(EXIT_FAILURE);
    }

    std::cout << std::endl;
    std::cout << "*** Reading meshfile info ***" << std::endl;

    /*
     * recall that inside mesh_info.dat the units are:
     * theta:  radians
     * radius: km
     * so for the radius, we need to convert to meters
     */
    std::cout << std::endl;
    std::ifstream source;
    source.open( filePath, std::ios_base::in);
    std::string line;
    while (std::getline(source, line) )
      {
        std::istringstream ss(line);
        // colVal is used to store the value of the column being read
        std::string colVal;
        // read first column
        ss >> colVal;

        if (colVal == "thL"){
          ss >> colVal; domainBounds_[0] = std::stod(colVal);
          std::cout << "thetaLeft (deg) = "
		    << std::setprecision(dblFmt)
		    << domainBounds_[0] << std::endl;
        }

        if (colVal == "thR"){
          ss >> colVal; domainBounds_[1] = std::stod(colVal);
          std::cout << "thetaRight (deg) = "
		    << std::setprecision(dblFmt)
		    << domainBounds_[1] << std::endl;
        }

        if (colVal == "rCmb"){
          // multiply by 1000 to convert from km to m
          ss >> colVal; domainBounds_[2] = std::stod(colVal)*thousand;
          std::cout << "minimum radius (km) = "
		    << std::setprecision(dblFmt)
		    << domainBounds_[2]/thousand << std::endl;
        }

        if (colVal == "rSurf"){
          // multiply by 1000 to convert from km to m
          ss >> colVal; domainBounds_[3] = std::stod(colVal)*thousand;
          std::cout << "surface radius (km) = "
		    << std::setprecision(dblFmt)
		    << domainBounds_[3]/thousand << std::endl;
        }

        if (colVal == "dth"){
          ss >> colVal; dth_ = std::stod(colVal);
          dthInv_ = one/dth_;
          std::cout << "dth [rad] = "
		    << std::setprecision(dblFmt)
		    << dth_ << std::endl;
        }

        if (colVal == "dr"){
          // multiply by 1000 to convert from km to m
          ss >> colVal; drr_ = std::stod(colVal)*thousand;
          drrInv_ = one/drr_;
          std::cout << "drr [km] = "
		    << std::setprecision(dblFmt)
		    << drr_/thousand << std::endl;
        }

        if (colVal == "numPtsVp"){
          ss >> colVal; numGptVp_ = std::stoi(colVal);
          std::cout << "numGptVp = " << numGptVp_ << std::endl;
        }
        if (colVal == "numPtsSp"){
          ss >> colVal; numGptSp_ = std::stoi(colVal);
          std::cout << "numGptSp = " << numGptSp_ << std::endl;
        }

        if (colVal == "nth"){
          ss >> colVal; numPtsAlongTh_ = std::stoi(colVal);
          std::cout << "nth = " << numPtsAlongTh_ << std::endl;
        }
        if (colVal == "nr"){
          ss >> colVal; numPtsAlongR_ = std::stoi(colVal);
          std::cout << "nr = " << numPtsAlongR_ << std::endl;
        }
      }//while
  }//readMeshInfo

};
#endif
