
#ifndef SHAXIPP_KOKKOS_TYPES_HPP_
#define SHAXIPP_KOKKOS_TYPES_HPP_

namespace kokkosapp{

struct commonTypes
{
  using scalar_type = double;

  using p_gs_t  = ParserGeneralSection<scalar_type>;
  using p_io_t  = ParserIoSection<scalar_type>;
  using p_mm_t  = ParserMaterialModel<scalar_type>;
  using p_ss_t  = ParserForcingSection<scalar_type>;
  using parser_type  = InputParser<p_gs_t, p_io_t, p_mm_t, p_ss_t>;

  using mesh_info_type = MeshInfo<scalar_type>;

  using exe_space = Kokkos::DefaultExecutionSpace;

  // jacobian is a sparse matrix
  using jacobian_ord_type = typename mesh_info_type::ordinal_type;
  using jacobian_d_type = KokkosSparse::CrsMatrix<scalar_type, jacobian_ord_type, exe_space>;

  // observer
  using observer_type = StateObserver<scalar_type>;

  // seismogram
  using seismogram_type  = Seismogram<scalar_type>;
};

struct rank1Types : commonTypes
{
  using typename commonTypes::scalar_type;
  using typename commonTypes::parser_type;
  using typename commonTypes::mesh_info_type;
  using typename commonTypes::jacobian_ord_type;
  using typename commonTypes::jacobian_d_type;
  using typename commonTypes::observer_type;
  using typename commonTypes::seismogram_type;

  // state is a rank-1 view
  using state_d_type = Kokkos::View<scalar_type*, Kokkos::DefaultExecutionSpace>;
  using state_h_type = typename state_d_type::host_mirror_type;

  // forcing type
  using forcing_type = RankOneForcing<scalar_type, state_d_type>;
};

struct rank2Types : commonTypes
{
  using typename commonTypes::scalar_type;
  using typename commonTypes::parser_type;
  using typename commonTypes::mesh_info_type;
  using typename commonTypes::jacobian_ord_type;
  using typename commonTypes::jacobian_d_type;
  using typename commonTypes::observer_type;
  using typename commonTypes::seismogram_type;

  // state is a rank-2 view
  using state_d_type = Kokkos::View<scalar_type**, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace>;
  using state_h_type = typename state_d_type::host_mirror_type;

  // forcing type
  using signal_type = Signal<scalar_type>;
  using signal_instances_h_type = Kokkos::View<signal_type*, Kokkos::HostSpace>;
  using forcing_type = RankTwoForcing<scalar_type, signal_instances_h_type>;
};

}// end namespace
#endif
