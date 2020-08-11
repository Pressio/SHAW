
#ifndef SHWAVEPP_KOKKOS_HPP_
#define SHWAVEPP_KOKKOS_HPP_

namespace kokkosapp{

template <
  typename sc_t,
  typename int_t,
  typename mesh_info_t,
  typename jacobian_d_t,
  typename exespace
  >
class ShWavePP
{
public:
  using scalar_type = sc_t;
  using klr = Kokkos::LayoutRight;
  using kll = Kokkos::LayoutLeft;

  using jacobian_ord_type = typename jacobian_d_t::ordinal_type;

  using coords_d_t	= Kokkos::View<sc_t*[2], klr, exespace>;
  using coords_h_t	= typename coords_d_t::host_mirror_type;

  using graph_vp_d_t	= Kokkos::View<int_t*[5], klr, exespace>;
  using graph_vp_h_t	= typename graph_vp_d_t::host_mirror_type;
  using graph_sp_d_t	= Kokkos::View<int_t*[3], klr, exespace>;
  using graph_sp_h_t	= typename graph_sp_d_t::host_mirror_type;

  using sten_coeff_d_t	= Kokkos::View<sc_t*[4], klr, exespace>;
  using sten_coeff_h_t	= typename sten_coeff_d_t::host_mirror_type;

  using rho_inv_d_t	= Kokkos::View<sc_t*, kll, exespace>;
  using rho_inv_h_t	= typename rho_inv_d_t::host_mirror_type;
  using rho_d_t		= Kokkos::View<sc_t*, kll, exespace>;
  using rho_h_t		= typename rho_d_t::host_mirror_type;
  using shmod_d_t	= Kokkos::View<sc_t*, kll, exespace>;
  using shmod_h_t	= typename shmod_d_t::host_mirror_type;
  using cot_d_t		= Kokkos::View<sc_t*, kll, exespace>;
  using cot_h_t		= typename cot_d_t::host_mirror_type;
  using labels_d_t	= Kokkos::View<int_t*, kll, exespace>;
  using labels_h_t	= typename labels_d_t::host_mirror_type;

  // static constexpr auto zero	= constants<sc_t>::zero();
  static constexpr auto one	= constants<sc_t>::one();
  static constexpr auto two	= constants<sc_t>::two();
  static constexpr auto three	= constants<sc_t>::three();
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
    this->allocateStorage();
  }

  void computeJacobiansWithoutMatProp(const MaterialModelBase<sc_t> & matModel)
  {
    std::cout << std::endl;
    std::cout << "*** Compute FOM Jacobians: without mat prop ***" << std::endl;

    sten_coeff_h_t coeffsVp_h("stenCoeffVp", numGptVp_);
    cot_h_t cotVp_h("cotVph", numGptVp_);
    cot_h_t cotSp_h("cotSph", numGptSp_);

    // read Vp graphs
    readFullMeshGraphFile<sc_t>(meshDir_, dofId::vp, graphVp_h_, coordsVp_h_, cotVp_h);
    readFullMeshCoeffFile<sc_t>(meshDir_, dofId::vp, coeffsVp_h);

    // 2. repeate for Sp
    readFullMeshGraphFile<sc_t>(meshDir_, dofId::sp, graphSp_h_, coordsSp_h_,
				cotSp_h, labelsSp_h_);

    fillVpJacobian(cotVp_h, coeffsVp_h);
    fillSpJacobian(cotSp_h);

    //after filling Jacobians, we can store material properties
    this->setMaterialProperties(matModel);

    printJacInfo();
  }

  void computeJacobiansWithMatProp(const MaterialModelBase<sc_t> & matModel)
  {
    std::cout << std::endl;
    std::cout << "*** Compute FOM Jacobians: with mat prop ***" << std::endl;

    sten_coeff_h_t coeffsVp_h("stenCoeffVp", numGptVp_);
    cot_h_t cotVp_h("cotVph", numGptVp_);
    cot_h_t cotSp_h("cotSph", numGptSp_);

    // read Vp graphs
    readFullMeshGraphFile<sc_t>(meshDir_, dofId::vp, graphVp_h_, coordsVp_h_, cotVp_h);
    readFullMeshCoeffFile<sc_t>(meshDir_, dofId::vp, coeffsVp_h);

    // 2. repeate for Sp
    readFullMeshGraphFile<sc_t>(meshDir_, dofId::sp, graphSp_h_, coordsSp_h_,
				cotSp_h, labelsSp_h_);

    // store material properties before filling jacs since we need them
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

  auto viewGraphHost(const dofId dof) const{
    switch(dof){
    case dofId::vp: return graphVp_h_; break;
    case dofId::sp: return graphSp_h_; break;
    default: throw std::runtime_error("Invalid dof");
    }
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

  sc_t getMinShearWaveVelocity() const{ return minMaxShearWaveVelocity_[0]; }
  sc_t getMaxShearWaveVelocity() const{ return minMaxShearWaveVelocity_[1]; }

  void writeCoordinatesToFile(dofId dof)
  {
    const std::string dofName = dofIdToString(dof);
    const std::string filePath = "coords_" + dofName + ".txt";
    const auto coords = this->viewCoordsHost(dof);
    std::ofstream file; file.open(filePath);
    for(auto i=0; i < coords.extent(0); i++){
      file << std::setprecision(dblFmt)
	   << coords(i,0) << " " << constants<sc_t>::one()/coords(i,1)
	   << std::endl;
    }
    file.close();
  }

private:
  void setMaterialProperties(const MaterialModelBase<sc_t> & matModel)
  {
    rhoInvVp_h_   = Kokkos::create_mirror_view(rhoInvVp_d_);
    Kokkos::deep_copy(rhoInvVp_h_, rhoInvVp_d_);
    auto shearModSp_h_ = Kokkos::create_mirror_view(shearModSp_d_);
    Kokkos::deep_copy(shearModSp_h_, shearModSp_d_);

    const auto gidsVp = this->viewGidListHost(dofId::vp);
    const auto gidsSp = this->viewGidListHost(dofId::sp);
    sc_t tmpRho, tmpVs;

    // set prop for the vp grid points
    for (int_t iPt=0; iPt < numGptVp_; ++iPt){
      const auto ptGID        = gidsVp(iPt);
      const auto thisPtTheta  = coordsVp_h_(ptGID, 0);
      const auto thisPtRadius = one/coordsVp_h_(ptGID, 1);
      matModel.computeAt(thisPtRadius, thisPtTheta, tmpRho, tmpVs);
      rhoInvVp_h_(iPt) = one/tmpRho;
    }

    // set prop for the stress grid points
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

    std::cout << "minMaxVs = "
	      << minMaxShearWaveVelocity_[0] << " "
	      << minMaxShearWaveVelocity_[1]
	      << std::endl;

    Kokkos::deep_copy(rhoInvVp_d_,   rhoInvVp_h_);
    Kokkos::deep_copy(shearModSp_d_, shearModSp_h_);
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

  auto countNNZ() const
  {
    const int_t numRows = numGptVp_;
    int_t result = 0;
    for (int_t iPt=0; iPt < numRows; ++iPt){
      const auto & ptGID      = graphVp_h_(iPt, 0);
      const auto & gid_west   = graphVp_h_(iPt, 1);
      const auto & gid_north  = graphVp_h_(iPt, 2);
      const auto & gid_east   = graphVp_h_(iPt, 3);
      const auto & gid_south  = graphVp_h_(iPt, 4);

      int_t thisRowNnz = 0;
      if (gid_north == gid_south)
	thisRowNnz+=1;
      else
	thisRowNnz+=2;

      if (gid_west == gid_east)
	thisRowNnz+=1;
      else
	thisRowNnz+=2;

      result+=thisRowNnz;
    }
    return result;
  }

  void fillVpJacobian(const cot_h_t cotVp_h,
		      const sten_coeff_h_t coeffsVp_h,
		      bool includeMatProp = false)
  {
    auto rhoInvVp_h = Kokkos::create_mirror_view(rhoInvVp_d_);

    const int_t numRows = numGptVp_;
    const int_t numCols = numGptSp_;

    // count nnz
    const auto nnz = countNNZ();

    // create data on device that we need to fill to create Jacobian
    typename jacobian_d_t::row_map_type::non_const_type ptr ("ptr", numRows+1);
    typename jacobian_d_t::index_type::non_const_type ind ("ind", nnz);
    typename jacobian_d_t::values_type val ("val", nnz);

    {
      // create host mirros of these
      typename jacobian_d_t::row_map_type::HostMirror ptr_h = Kokkos::create_mirror_view (ptr);
      typename jacobian_d_t::index_type::HostMirror   ind_h = Kokkos::create_mirror_view (ind);
      typename jacobian_d_t::values_type::HostMirror  val_h = Kokkos::create_mirror_view (val);

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

	const auto rhoInvOpt	= includeMatProp ? rhoInvVp_h[iPt] : one;

	const auto & c0	    = coeffsVp_h(iPt, 0);
	const auto & c1	    = coeffsVp_h(iPt, 1);
	const auto & c2	    = coeffsVp_h(iPt, 2);
	const auto & c3	    = coeffsVp_h(iPt, 3);

	const auto c_west  = (-rInv*dthInv_ + rInv*cotVp_h[ptGID])*c0*rhoInvOpt;
	const auto c_north = (drrInv_	    + three*oneHalf*rInv )*c1*rhoInvOpt;
	const auto c_east  = (rInv*dthInv_  + rInv*cotVp_h[ptGID])*c2*rhoInvOpt;
	const auto c_south = (-drrInv_	    + three*oneHalf*rInv )*c3*rhoInvOpt;

	int_t shift = 0;
	if (gid_north == gid_south){
	  shift+=1;
	  ind_h[++k] = gid_north;
	  val_h[k] = c_north+c_south;
	}
	else{
	  shift+=2;
	  ind_h[++k] = gid_north;
	  val_h[k] = c_north;

	  ind_h[++k] = gid_south;
	  val_h[k] = c_south;
	}

	if (gid_west == gid_east){
	  shift+=1;
	  ind_h[++k] = gid_west;
	  val_h[k] = c_west+c_east;
	}
	else{
	  shift+=2;
	  ind_h[++k] = gid_west;
	  val_h[k] = c_west;

	  ind_h[++k] = gid_east;
	  val_h[k] = c_east;
	}
	ptr_h[iPt+1] = ptr_h[iPt] + shift;
      }
      Kokkos::deep_copy(ptr, ptr_h);
      Kokkos::deep_copy(ind, ind_h);
      Kokkos::deep_copy(val, val_h);
    }

    jacobian_d_t J("JacVp", numRows, numCols, nnz, val, ptr, ind);
    JacVp_d_ = J;
  }

  void fillSpJacobian(const cot_h_t cotSp_h, bool includeMatProp = false)
  {
    auto shearModSp_h = Kokkos::create_mirror_view(shearModSp_d_);

    const int_t nonZerosPerRowJ_ = 2;
    const int_t numRows = numGptSp_;
    const int_t numCols = numGptVp_;
    const int_t numEnt  = numRows * nonZerosPerRowJ_;

    // create data on device that we need to fill to create Jacobian
    typename jacobian_d_t::row_map_type::non_const_type ptr ("ptr", numRows+1);
    typename jacobian_d_t::index_type::non_const_type ind ("ind", numEnt);
    typename jacobian_d_t::values_type val ("val", numEnt);

    {
      // create host mirros of these
      typename jacobian_d_t::row_map_type::HostMirror ptr_h = Kokkos::create_mirror_view (ptr);
      typename jacobian_d_t::index_type::HostMirror   ind_h = Kokkos::create_mirror_view (ind);
      typename jacobian_d_t::values_type::HostMirror  val_h = Kokkos::create_mirror_view (val);

      // first, fill in how many elements per row
      ptr_h[0] = 0;
      for (int_t iRow = 0; iRow < numRows; ++iRow) {
    	ptr_h[iRow+1] = ptr_h[iRow] + nonZerosPerRowJ_;
      }
      Kokkos::deep_copy(ptr, ptr_h);

      // set the column index for each non-zero entry of the Jacobian
      int_t k = -1;
      for (int_t iPt=0; iPt < numRows; ++iPt)
      {
	const auto & ptGID       = graphSp_h_(iPt, 0);
	const auto & thisPtTheta = coordsSp_h_(ptGID, 0);
	const auto & rInv        = coordsSp_h_(ptGID, 1);
	const auto & myLabel	 = labelsSp_h_(iPt);

	const auto shearModOpt	 = includeMatProp ? shearModSp_h[iPt] : one;

	if( myLabel==1){ //srp
	  const auto gid_north     = graphSp_h_(iPt, 1);
	  const auto gid_south     = graphSp_h_(iPt, 2);
	  const auto c_north	   = (drrInv_ - rInv*oneHalf)*shearModOpt;
	  const auto c_south	   = (-drrInv_ - rInv*oneHalf)*shearModOpt;
	  ind_h[++k] = gid_north;
	  val_h[k]   = c_north;
	  ind_h[++k] = gid_south;
	  val_h[k]   = c_south;
	}

	if(myLabel==2){ //stp
    	  const auto gid_west = graphSp_h_(iPt, 1);
    	  const auto gid_east = graphSp_h_(iPt, 2);
    	  const auto c_west   = (-dthInv_ - oneHalf*cotSp_h[ptGID])*rInv*shearModOpt;
    	  const auto c_east   = (dthInv_ - oneHalf*cotSp_h[ptGID])*rInv*shearModOpt;
    	  ind_h[++k] = gid_west;
    	  val_h[k]   = c_west;
    	  ind_h[++k] = gid_east;
    	  val_h[k]   = c_east;
	}
      }

      Kokkos::deep_copy(ind, ind_h);
      Kokkos::deep_copy(val, val_h);
    }

    jacobian_d_t J("JacSp", numRows, numCols, numEnt, val, ptr, ind);
    JacSp_d_ = J;
  }

  void allocateStorage()
  {
    // we don't need to pre-allocate the jacobians because those
    // are created during assemble

    // for Vp
    Kokkos::resize(coordsVp_h_,	numGptVp_);
    Kokkos::resize(rhoInvVp_d_, numGptVp_);
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
  sc_t dthInv_{};
  sc_t drrInv_{};

  // bool to tell if the jacobians include the material prop
  // or if the mat properties are factored out
  bool factorMatPropFromJac_ = {};

  // min and max value of the shear wave velocity
  std::array<sc_t, 2> minMaxShearWaveVelocity_ = {std::numeric_limits<sc_t>::max(),
						  std::numeric_limits<sc_t>::min()};

  //**************************
  //***** members for Vp *****
  //**************************
  // numGptVp = number of points for the velocity (Vp)
  int_t numGptVp_		= {};

  // coords for velocity point (we store theta and r) - host only
  // theta in col[0], 1/r in col[1]
  coords_h_t coordsVp_h_ = {};

  // vector containing 1/density for velocity points
  rho_inv_d_t rhoInvVp_d_	= {};
  rho_inv_d_t rhoInvVp_h_	= {};

  // graph for vp
  graph_vp_h_t  graphVp_h_ = {};

  // jacobian matrix for Vp
  jacobian_d_t JacVp_d_	= {};

  //**************************
  //***** members for Sp *****
  //**************************
  // numGptSp_ = number of stress points
  int_t numGptSp_ = {};

  // coords for velocity point (we store theta and r) - host only
  // theta in col[0], 1/r in col[1]
  coords_h_t coordsSp_h_ = {};

  labels_h_t labelsSp_h_ = {};

  // density and shear modulus for sp points
  shmod_d_t shearModSp_d_	= {};

  // graph for sp
  graph_sp_h_t  graphSp_h_ = {};

  // jacobian matrix for sp
  jacobian_d_t JacSp_d_ = {};

};//end class

}//end namespace kokkosapp
#endif
