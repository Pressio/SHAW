
#include "../shared/all.hpp"
#include "common_types.hpp"
#include "shwavepp.hpp"
#include "./fom/run_fom.hpp"
#include "./fom/fom_problem_rank_one_forcing.hpp"
// #include "./fom/fom_problem_rank_two_forcing.hpp"

template<typename scalar_t>
struct MyCustomMaterialModel final : public MaterialModelBase<scalar_t>
{

  template<typename mesh_info_t>
  MyCustomMaterialModel(const mesh_info_t & meshInfo)
  {
    // constuct how you want
  }

  void computeAt(const scalar_t & radiusFromCenterMeters,
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

    // create parser for yaml file
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

    // if(parser.enableMultiForcing())
    // {
    //   // using prob_t = kokkosapp::FomProblemRankTwoForcing;
    //   // prob_t problem(parser, meshInfo, materialModel);
    //   // problem.execute();
    // }
    // else{
      using prob_t = kokkosapp::FomProblemRankOneForcing<kokkosapp::rank1Types>;
      prob_t problem(parser, meshInfo, materialModel);
      problem.execute();
    // }
  }
  Kokkos::finalize();

  return 0;
}
