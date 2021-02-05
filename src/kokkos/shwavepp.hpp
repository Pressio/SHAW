
#ifndef SHWAVEPP_KOKKOS_HPP_
#define SHWAVEPP_KOKKOS_HPP_

namespace kokkosapp{

template <
  typename scalar_type,
  typename int_t,
  typename mesh_info_t,
  typename jacobian_d_t,
  typename exespace
  >
class ShWavePP
{

public:
  using klr = Kokkos::LayoutRight;
  using kll = Kokkos::LayoutLeft;

  using jacobian_ord_type = typename jacobian_d_t::ordinal_type;

  using coords_h_t	= Kokkos::View<scalar_type*[2], klr, Kokkos::HostSpace>;

  using graph_vp_d_t	= Kokkos::View<int_t*[5], klr, exespace>;
  using graph_vp_h_t	= typename graph_vp_d_t::host_mirror_type;
  using graph_sp_d_t	= Kokkos::View<int_t*[3], klr, exespace>;
  using graph_sp_h_t	= typename graph_sp_d_t::host_mirror_type;

  using sten_coeff_h_t	= Kokkos::View<scalar_type*[4], klr, Kokkos::HostSpace>;

  using rho_inv_d_t	= Kokkos::View<scalar_type*, exespace>;
  using rho_inv_h_t	= typename rho_inv_d_t::host_mirror_type;
  using rho_d_t		= Kokkos::View<scalar_type*, exespace>;
  using rho_h_t		= typename rho_d_t::host_mirror_type;
  using shmod_d_t	= Kokkos::View<scalar_type*, exespace>;
  using shmod_h_t	= typename shmod_d_t::host_mirror_type;
  using cot_d_t		= Kokkos::View<scalar_type*, exespace>;
  using cot_h_t		= typename cot_d_t::host_mirror_type;
  using labels_d_t	= Kokkos::View<int_t*, exespace>;
  using labels_h_t	= typename labels_d_t::host_mirror_type;

  static constexpr auto one	= constants<scalar_type>::one();
  static constexpr auto two	= constants<scalar_type>::two();
  static constexpr auto three	= constants<scalar_type>::three();
  static constexpr auto oneHalf = one/two;

public:
  ShWavePP() = delete;

  ShWavePP(const mesh_info_t & meshInfo)
    : meshDir_{meshInfo.getMeshDir()},
      dthInv_{meshInfo.getAngularSpacingInverse()},
      drrInv_{meshInfo.getRadialSpacingInverse()},
      numGptVp_{meshInfo.getNumVpPts()},
      numGptSp_{meshInfo.getNumSpPts()}
  {
    this->allocateMembers();
  }

  void computeJacobians(const MaterialModelBase<scalar_type> & matModel)
  {
    std::cout << std::endl;
    std::cout << "*** Compute FOM Jacobian matrices ***" << std::endl;

    sten_coeff_h_t coeffsVp_h("stenCoeffVp", numGptVp_);
    cot_h_t cotVp_h("cotVph", numGptVp_);
    cot_h_t cotSp_h("cotSph", numGptSp_);

    // read Vp graph
    readFullMeshGraphFile<scalar_type>(meshDir_, dofId::vp, graphVp_h_, coordsVp_h_, cotVp_h);
    readFullMeshCoeffFile<scalar_type>(meshDir_, dofId::vp, coeffsVp_h);

    // read Sp graph
    readFullMeshGraphFile<scalar_type>(meshDir_, dofId::sp, graphSp_h_, coordsSp_h_,
				       cotSp_h, labelsSp_h_);

    // store material properties since are needed to fill jacobians
    this->setMaterialProperties(matModel);

    fillVpJacobian(cotVp_h, coeffsVp_h, true);
    fillSpJacobian(cotSp_h, true);

    printJacInfo();
  }

  const auto viewGidListHost(const dofId dof) const{
    switch(dof){
    case dofId::vp:
      return Kokkos::subview(graphVp_h_, Kokkos::ALL(), 0);
      break;
    case dofId::sp:
      return Kokkos::subview(graphSp_h_, Kokkos::ALL(), 0);
      break;
    default:
      throw std::runtime_error("Invalid dof");
    }
  }

  auto viewVelocityGraphHost() const{
    return graphVp_h_;
  }

  auto viewStressGraphHost() const{
    return graphSp_h_;
  }

  auto viewCoordsHost(const dofId dof) const{
    switch(dof){
    case dofId::vp: return coordsVp_h_; break;
    case dofId::sp: return coordsSp_h_; break;
    default: throw std::runtime_error("Invalid dof");
    }
  }

  auto viewLabelsHost(const dofId dof) const{
    switch(dof){
    case dofId::sp: return labelsSp_h_; break;
    default: throw std::runtime_error("Invalid dof");
    }
  }

  auto viewJacobianDevice(const dofId dof) const{
    switch(dof){
    case dofId::vp: return JacVp_d_; break;
    case dofId::sp: return JacSp_d_; break;
    default: throw std::runtime_error("Invalid dof");
    }
  }

  auto viewInvDensityDevice(const dofId dof) const{
    switch(dof){
    case dofId::vp: return rhoInvVp_d_; break;
    default: throw std::runtime_error("Invalid dof");
    }
  }

  auto viewInvDensityHost(const dofId dof) const{
    switch(dof){
    case dofId::vp: return rhoInvVp_h_; break;
    default: throw std::runtime_error("Invalid dof");
    }
  }

  auto viewShearModulusDevice(const dofId dof) const{
    switch(dof){
    case dofId::sp: return shearModSp_d_; break;
    default: throw std::runtime_error("Invalid dof");
    }
  }

  auto getJacobianNNZ(const dofId dof) const{
    switch(dof){
    case dofId::vp: return JacVp_d_.nnz(); break;
    case dofId::sp: return JacSp_d_.nnz(); break;
    default: throw std::runtime_error("Invalid dof");
    }
  }

  scalar_type getMinShearWaveVelocity() const{ return minMaxShearWaveVelocity_[0]; }
  scalar_type getMaxShearWaveVelocity() const{ return minMaxShearWaveVelocity_[1]; }

  void writeCoordinatesToFile(dofId dof)
  {
    const std::string dofName = dofIdToString(dof);
    const std::string filePath = "coords_" + dofName + ".txt";
    const auto coords = this->viewCoordsHost(dof);
    std::ofstream file; file.open(filePath);
    for(auto i=0; i < coords.extent(0); i++){
      file << std::setprecision(dblFmt)
	   << coords(i,0) << " " << constants<scalar_type>::one()/coords(i,1)
	   << std::endl;
    }
    file.close();
  }

private:
  void setMaterialProperties(const MaterialModelBase<scalar_type> & matModel)
  {
    rhoInvVp_h_   = Kokkos::create_mirror_view(rhoInvVp_d_);
    shearModSp_h_ = Kokkos::create_mirror_view(shearModSp_d_);
    KokkosBlas::fill(rhoInvVp_h_,   static_cast<scalar_type>(0));
    KokkosBlas::fill(shearModSp_h_, static_cast<scalar_type>(0));

    const auto gidsVp = this->viewGidListHost(dofId::vp);
    const auto gidsSp = this->viewGidListHost(dofId::sp);
    scalar_type tmpRho, tmpVs;

    // set host properties for the vp dofs
    for (int_t iPt=0; iPt < numGptVp_; ++iPt){
      const auto ptGID        = gidsVp(iPt);
      const auto thisPtTheta  = coordsVp_h_(ptGID, 0);
      const auto thisPtRadius = one/coordsVp_h_(ptGID, 1);
      matModel.computeAt(thisPtRadius, thisPtTheta, tmpRho, tmpVs);
      rhoInvVp_h_(iPt) = one/tmpRho;
    }
    Kokkos::deep_copy(rhoInvVp_d_, rhoInvVp_h_);

    // set host properties for the stress dofs
    for (int_t iPt=0; iPt < numGptSp_; ++iPt){
      const auto ptGID	      = gidsSp(iPt);
      const auto thisPtTheta  = coordsSp_h_(ptGID, 0);
      const auto thisPtRadius = one/coordsSp_h_(ptGID, 1);
      matModel.computeAt(thisPtRadius, thisPtTheta, tmpRho, tmpVs);
      shearModSp_h_(iPt) = tmpRho * tmpVs * tmpVs;

      if (tmpVs != 0.){
    	minMaxShearWaveVelocity_[0] = std::min(minMaxShearWaveVelocity_[0], tmpVs);
      }
      minMaxShearWaveVelocity_[1] = std::max(minMaxShearWaveVelocity_[1], tmpVs);
    }
    Kokkos::deep_copy(shearModSp_d_, shearModSp_h_);

    std::cout << "minMaxVs = "
    	      << minMaxShearWaveVelocity_[0] << " "
    	      << minMaxShearWaveVelocity_[1]
    	      << std::endl;
  }

  auto countVpJacNNZ() const
  {
    int_t result = 0;
    for (int_t iPt=0; iPt < numGptVp_; ++iPt)
    {
      const auto & ptGID      = graphVp_h_(iPt, 0);
      const auto & gid_west   = graphVp_h_(iPt, 1);
      const auto & gid_north  = graphVp_h_(iPt, 2);
      const auto & gid_east   = graphVp_h_(iPt, 3);
      const auto & gid_south  = graphVp_h_(iPt, 4);

      int_t thisRowNnz = 0;
      if (gid_north == gid_south){
	thisRowNnz+=1;
      }
      else{
	thisRowNnz+=2;
      }

      if (gid_west == gid_east){
	thisRowNnz+=1;
      }
      else{
	thisRowNnz+=2;
      }

      result+=thisRowNnz;
    }
    return result;
  }

  void fillVpJacobian(const cot_h_t cotVp_h,
		      const sten_coeff_h_t coeffsVp_h,
		      bool includeMatProp = false)
  {
    const int_t numRows = numGptVp_;
    const int_t numCols = numGptSp_;

    // count nnz
    const auto nnz = countVpJacNNZ();

    // // create data on device that we need to fill to create Jacobian
    // typename jacobian_d_t::values_type val ("val", nnz);
    // typename jacobian_d_t::row_map_type::non_const_type ptr ("ptr", numRows+1);
    // typename jacobian_d_t::index_type::non_const_type ind ("ind", nnz);
    // // create host mirros of these
    // auto val_h = Kokkos::create_mirror_view (val);
    // auto ptr_h = Kokkos::create_mirror_view (ptr);
    // auto ind_h = Kokkos::create_mirror_view (ind);

    Kokkos::View<scalar_type*, Kokkos::HostSpace> val_h("valVp", nnz);
    Kokkos::View<int_t*, Kokkos::HostSpace> ptr_h("ptrVp", numRows+1);
    Kokkos::View<int_t*, Kokkos::HostSpace> ind_h("ptrVp", nnz);

    ptr_h[0] = 0;
    // set the column index for each non-zero entry of the Jacobian
    int_t k = -1;
    for (int_t iPt=0; iPt < numRows; ++iPt)
      {
	const auto & ptGID	 = graphVp_h_(iPt, 0);
	const auto & gid_west	 = graphVp_h_(iPt, 1);
	const auto & gid_north   = graphVp_h_(iPt, 2);
	const auto & gid_east    = graphVp_h_(iPt, 3);
	const auto & gid_south   = graphVp_h_(iPt, 4);
	const auto & thisPtTheta = coordsVp_h_(ptGID, 0);
	const auto & rInv	 = coordsVp_h_(ptGID, 1);
	const auto rhoInv	 = rhoInvVp_h_[iPt];

	const auto & c0	    = coeffsVp_h(iPt, 0);
	const auto & c1	    = coeffsVp_h(iPt, 1);
	const auto & c2	    = coeffsVp_h(iPt, 2);
	const auto & c3	    = coeffsVp_h(iPt, 3);

	// this is here because it is used for a test to check correctness of entries
	const auto c_west  = (-rInv*dthInv_ + rInv*cotVp_h[ptGID])*c0*rhoInv;
	const auto c_north = (drrInv_	    + three*oneHalf*rInv )*c1*rhoInv;
	const auto c_east  = (rInv*dthInv_  + rInv*cotVp_h[ptGID])*c2*rhoInv;
	const auto c_south = (-drrInv_	    + three*oneHalf*rInv )*c3*rhoInv;

	int_t shift = 0;
	if (gid_north == gid_south){
	  shift+=1;
	  ind_h[++k] = gid_north;
	  val_h(k) = c_north+c_south;
	}
	else{
	  shift+=2;
	  ind_h[++k] = gid_north;
	  val_h(k) = c_north;

	  ind_h[++k] = gid_south;
	  val_h(k) = c_south;
	}

	if (gid_west == gid_east){
	  shift+=1;
	  ind_h[++k] = gid_west;
	  val_h(k) = c_west+c_east;
	}
	else{
	  shift+=2;
	  ind_h[++k] = gid_west;
	  val_h(k) = c_west;

	  ind_h[++k] = gid_east;
	  val_h(k) = c_east;
	}
	ptr_h[iPt+1] = ptr_h[iPt] + shift;
      }

    jacobian_d_t J("JacVp", numRows, numCols, nnz,
		   val_h.data(), ptr_h.data(), ind_h.data());
    JacVp_d_ = J;

    // auto h_e = Kokkos::create_mirror_view(JacVp_d_.graph.entries);
    // Kokkos::deep_copy(h_e, JacVp_d_.graph.entries);
    // for (auto i=0; i<h_e.extent(0); ++i)
    //   std::cout << h_e(i) << std::endl;

    // auto h_r = Kokkos::create_mirror_view(JacVp_d_.graph.row_map);
    // Kokkos::deep_copy(h_r, JacVp_d_.graph.row_map);
    // for (auto i=0; i<h_r.extent(0); ++i)
    //   std::cout << h_r(i) << std::endl;

    // auto h_v = Kokkos::create_mirror_view(JacVp_d_.values);
    // Kokkos::deep_copy(h_v, JacVp_d_.values);
    // for (auto i=0; i<h_v.extent(0); ++i)
    //   std::cout << h_v(i) << std::endl;
  }

  void fillSpJacobian(const cot_h_t cotSp_h, bool includeMatProp = false)
  {
    //auto shearModSp_h = Kokkos::create_mirror_view(shearModSp_d_);

    const int_t nonZerosPerRowJ_ = 2;
    const int_t numRows = numGptSp_;
    const int_t numCols = numGptVp_;
    const int_t numEnt  = numRows * nonZerosPerRowJ_;

    // // create data on device that we need to fill to create Jacobian
    // typename jacobian_d_t::row_map_type::non_const_type ptr ("ptr", numRows+1);
    // typename jacobian_d_t::index_type::non_const_type ind ("ind", numEnt);
    // typename jacobian_d_t::values_type val ("val", numEnt);
    // // create host mirros of these
    // typename jacobian_d_t::row_map_type::HostMirror ptr_h = Kokkos::create_mirror_view (ptr);
    // typename jacobian_d_t::index_type::HostMirror   ind_h = Kokkos::create_mirror_view (ind);
    // typename jacobian_d_t::values_type::HostMirror  val_h = Kokkos::create_mirror_view (val);

    Kokkos::View<scalar_type*, Kokkos::HostSpace> val_h("valSp", numEnt);
    Kokkos::View<int_t*, Kokkos::HostSpace> ptr_h("ptrSp", numRows+1);
    Kokkos::View<int_t*, Kokkos::HostSpace> ind_h("indSp", numEnt);

    // first, fill in how many elements per row
    ptr_h[0] = 0;
    for (int_t iRow = 0; iRow < numRows; ++iRow) {
      ptr_h[iRow+1] = ptr_h[iRow] + nonZerosPerRowJ_;
    }

    // set the column index for each non-zero entry of the Jacobian
    int_t k = -1;
    for (int_t iPt=0; iPt < numRows; ++iPt)
      {
	const auto & ptGID       = graphSp_h_(iPt, 0);
	const auto & thisPtTheta = coordsSp_h_(ptGID, 0);
	const auto & rInv        = coordsSp_h_(ptGID, 1);
	const auto & myLabel	 = labelsSp_h_(iPt);
	const auto shearMod	 = shearModSp_h_(iPt);

	if( myLabel==1){ //srp
	  const auto gid_north     = graphSp_h_(iPt, 1);
	  const auto gid_south     = graphSp_h_(iPt, 2);
	  const auto c_north	   = (drrInv_ - rInv*oneHalf)*shearMod;
	  const auto c_south	   = (-drrInv_ - rInv*oneHalf)*shearMod;
	  ind_h[++k] = gid_north;
	  val_h[k]   = c_north;
	  ind_h[++k] = gid_south;
	  val_h[k]   = c_south;
	}

	if(myLabel==2){ //stp
    	  const auto gid_west = graphSp_h_(iPt, 1);
    	  const auto gid_east = graphSp_h_(iPt, 2);
    	  const auto c_west   = (-dthInv_ - oneHalf*cotSp_h[ptGID])*rInv*shearMod;
    	  const auto c_east   = (dthInv_ - oneHalf*cotSp_h[ptGID])*rInv*shearMod;
    	  ind_h[++k] = gid_west;
    	  val_h[k]   = c_west;
    	  ind_h[++k] = gid_east;
    	  val_h[k]   = c_east;
	}
      }

    jacobian_d_t J2("JacSp", numRows, numCols, numEnt,
		    val_h.data(), ptr_h.data(), ind_h.data());
    JacSp_d_ = J2;
  }

  void printJacInfo() const
  {
    std::cout << "jacVp: "
	      << " nnz = " << getJacobianNNZ(dofId::vp)
	      << " nrows = " << JacVp_d_.numRows()
	      << " ncols = " << JacVp_d_.numCols() << std::endl;

    std::cout << "jacSp: "
	      << " nnz = " << getJacobianNNZ(dofId::sp)
	      << " nrows = " << JacSp_d_.numRows()
	      << " ncols = " << JacSp_d_.numCols() << std::endl;
  }

  void allocateMembers()
  {
    // we don't need to allocate the jacobians because they
    // are allocated directly during assemble

    // for Vp
    Kokkos::resize(coordsVp_h_,	numGptVp_);
    Kokkos::resize(rhoInvVp_d_, numGptVp_);
    Kokkos::resize(rhoInvVp_h_, numGptVp_);
    Kokkos::resize(graphVp_h_,  numGptVp_);

    // for sp
    Kokkos::resize(coordsSp_h_,  numGptSp_);
    Kokkos::resize(labelsSp_h_,  numGptSp_);
    Kokkos::resize(shearModSp_d_,numGptSp_);
    Kokkos::resize(graphSp_h_,   numGptSp_);
  }

private:
  std::string meshDir_ = {};

  // inverse spacing in theta (rad) and r (m) direction
  scalar_type dthInv_{};
  scalar_type drrInv_{};

  // min and max value of the shear wave velocity
  std::array<scalar_type, 2> minMaxShearWaveVelocity_ =
    {std::numeric_limits<scalar_type>::max(),
     std::numeric_limits<scalar_type>::min()};

  //**************************
  //***** members for Vp *****
  //**************************
  // numGptVp = number of points for the velocity (Vp)
  int_t numGptVp_ = {};

  // coords for velocity point (we store theta and r) - host only
  // theta in col[0], 1/r in col[1]
  coords_h_t coordsVp_h_ = {};

  // array containing 1/density at each velocity point
  rho_inv_d_t rhoInvVp_d_ = {};
  rho_inv_h_t rhoInvVp_h_ = {};

  // graph for vp
  graph_vp_h_t  graphVp_h_ = {};

  // jacobian matrix for Vp
  jacobian_d_t JacVp_d_ = {};

  //**************************
  //***** members for Sp *****
  //**************************
  // numGptSp_ = number of stress points
  int_t numGptSp_ = {};

  // coords for velocity point (we store theta and r) - host only
  // theta in col[0], 1/r in col[1]
  coords_h_t coordsSp_h_ = {};

  // array storing labels to differentiate the stress dof,
  // i.e. sigma_r,phi from sigma_theta,phi
  labels_h_t labelsSp_h_ = {};

  // shear modulus for sp points
  shmod_d_t shearModSp_d_ = {};
  shmod_h_t shearModSp_h_ = {};

  // graph for sp
  graph_sp_h_t  graphSp_h_ = {};

  // jacobian matrix for sp
  jacobian_d_t JacSp_d_ = {};

};

}//end namespace kokkosapp
#endif
