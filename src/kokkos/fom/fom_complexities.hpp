
#ifndef FOM_COMPLEXITIES_HPP_
#define FOM_COMPLEXITIES_HPP_

namespace kokkosapp{

/*
  rank-1 forcing
*/
template <typename sc_t, typename state_d_t, typename app_t, typename forcing_t>
typename std::enable_if< is_kokkos_1dview<state_d_t>::value >::type
complexityFom(const state_d_t xVp,
	      const state_d_t xSp,
	      const app_t & fomObj,
	      const forcing_t & forcingObj,
	      double & memCostMB,
	      double & flopsCost)
{
  using ord_t = typename app_t::jacobian_ord_type;
  using comp_t = Complexity<sc_t>;

  const auto nVp = xVp.extent(0);
  const auto nSp = xSp.extent(0);

  std::array<double, 5> memMB = {};
  std::array<double, 5> flops = {};

  // spmv: xVp = xVp + dt * Jvp * xSp
  const auto nnz_j_vp = fomObj.getJacobianNNZ(dofId::vp);
  comp_t::template spmv<ord_t>(nnz_j_vp, nVp, memMB[1], flops[1]);

  // mult: xVp = xVp + dt * rhoInv * f   (note that beta=1)
  comp_t::mult_beta_one(nVp, memMB[2], flops[2]);

  // spmv: xSp = xSp + dt * Jsp * xVp
  const auto nnz_j_sp = fomObj.getJacobianNNZ(dofId::sp);
  comp_t::template spmv<ord_t>(nnz_j_sp, nSp, memMB[3], flops[3]);

  // accumulate
  memCostMB   = std::accumulate(memMB.begin(),   memMB.end(),   0.);
  flopsCost = std::accumulate(flops.begin(), flops.end(), 0.);
}


/*
   rank-2 forcing
*/
template <typename sc_t, typename state_d_t, typename app_t, typename forcing_t>
typename std::enable_if< is_kokkos_2dview<state_d_t>::value >::type
complexityFom(const state_d_t xVp,
	      const state_d_t xSp,
	      const app_t & fomObj,
	      const forcing_t & forcingObj,
	      double & memCostMB,
	      double & flopsCost)
{
  using ord_t = typename app_t::jacobian_ord_type;
  using comp_t = Complexity<sc_t>;

  assert( xVp.extent(1) == xSp.extent(1) );
  const auto nVp    = xVp.extent(0);
  const auto nSp    = xSp.extent(0);
  const auto fSize = xVp.extent(1);

  std::array<double, 5> memMB = {};
  std::array<double, 5> flops = {};

  // //account for forcing complexity
  // forcingObj.complexityOfEvaluateMethod(memMB[0], flops[0]);

  // xVp = xVp + dt * Jvp * xSp
  const auto nnz_j_vp = fomObj.getJacobianNNZ(dofId::vp);
  comp_t::template spmm<ord_t>(nnz_j_vp, nVp, fSize, memMB[1], flops[1]);

  memMB[2] = 4.*fSize*sizeof(sc_t);
  flops[2] = 3.*fSize;

  // spmm: xSp = xSp + dt * Jsp * xVp
  const auto nnz_j_sp = fomObj.getJacobianNNZ(dofId::sp);
  comp_t::template spmm<ord_t>(nnz_j_sp, nSp, fSize, memMB[3], flops[3]);

  // accumulate
  memCostMB = std::accumulate(memMB.begin(), memMB.end(), 0.);
  flopsCost = std::accumulate(flops.begin(), flops.end(), 0.);
}

}//end namespace kokkosapp
#endif
