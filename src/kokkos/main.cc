
// include order below matters

// common things
#include "../shared/all.hpp"
#include "types.hpp"
#include "shwavepp.hpp"
// fom
#include "fom_run.hpp"
#include "fom_problem_rank_one.hpp"
#include "fom_problem_rank_two.hpp"
// rom
#include "rom_load_basis.hpp"
#include "rom_compute_jacobians.hpp"
#include "rom_run_rank_one.hpp"
#include "rom_problem_rank_one.hpp"

template<typename scalar_t>
struct MyCustomMaterialModel final
  : public MaterialModelBase<scalar_t>
{

  template<typename mesh_info_t>
  MyCustomMaterialModel(const mesh_info_t & meshInfo){
    // as needed
  }

  void computeAt(const scalar_t & radiusFromCenterInMeters,
		 const scalar_t & angleRadians,
		 scalar_t & density,
		 scalar_t & vs) const final
  {
    // evaluate:
    // - density in units of kg/m^3
    // - shear velocity in units of [m/s]
    // at target grid point with coordinates (radiusFromCenterMeters, angleRadians)
    //...
  }
};

int main(int argc, char *argv[])
{
  Kokkos::initialize (argc, argv);
  {
    using scalar_t    = kokkosapp::commonTypes::scalar_type;
    using parser_t    = kokkosapp::commonTypes::parser_type;
    using mesh_info_t = kokkosapp::commonTypes::mesh_info_type;

    // create parser for input file
    parser_t parser(argc, argv);

    // create object with mesh info
    mesh_info_t meshInfo(parser.getMeshDir());

    // create material model object
    // if parser has custom model, use custom class above
    std::shared_ptr<MaterialModelBase<scalar_t>> materialModel = {};
    if (parser.getMaterialModelKind() == materialModelKind::custom){
      materialModel = std::make_shared<MyCustomMaterialModel<scalar_t>>(meshInfo);
    }
    else{
      materialModel = createMaterialModel<scalar_t>(parser, meshInfo);
    }

    // branch out for FOM or ROM
    if (parser.enableRom())
    {
      // if here, we want to run ROM
      using prob_t = kokkosapp::RomProblemRankOneForcing<kokkosapp::rank1TypesRom>;
      prob_t problem(parser, meshInfo, *materialModel);
      problem();

    }
    else
    {
      // if here, we want to run FOM

      if(parser.rank2Enabled())
	{
	  using prob_t = kokkosapp::FomProblemRankTwoForcing<kokkosapp::rank2Types>;
	  prob_t problem(parser, meshInfo, *materialModel);
	  problem();
	}
      else{
	using prob_t = kokkosapp::FomProblemRankOneForcing<kokkosapp::rank1Types>;
	prob_t problem(parser, meshInfo, *materialModel);
	problem();
      }
    }
  }
  Kokkos::finalize();

  return 0;
}
