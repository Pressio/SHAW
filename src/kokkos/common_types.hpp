
#ifndef SHAXIPP_KOKKOS_TYPES_HPP_
#define SHAXIPP_KOKKOS_TYPES_HPP_

namespace kokkosapp{

struct commonTypes
{
  using scalar_type	= double;
  using scalar_t	= scalar_type;
  using sc_t		= scalar_t;
  // type to use for all indexing, has to be large enough
  using int_t		= int64_t;

  using pgs_mx_t	= ParserGeneralSection<sc_t, int_t>;
  using pio_mx_t	= ParserIoSection<sc_t, int_t>;
  using pmm_mx_t	= ParserMaterialModel<sc_t, int_t>;
  using pss_mx_t	= ParserForcingSection<sc_t, int_t>;
  using prs_mx_t	= ParserRomSection<sc_t, int_t>;
  using psas_mx_t	= ParserSamplingSection<sc_t, int_t>;
  using parser_t	= InputParser<pgs_mx_t, pio_mx_t, pmm_mx_t, pss_mx_t, prs_mx_t, psas_mx_t>;

  // aliases for layouts and exe space
  using klr		= Kokkos::LayoutRight;
  using kll		= Kokkos::LayoutLeft;
  using exe_space	= Kokkos::DefaultExecutionSpace;
};

}// end namespace eigenapp
#endif
