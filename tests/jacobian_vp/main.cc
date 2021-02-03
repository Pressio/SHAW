
#include "./fom/fom_problem_rank_one_forcing.hpp"

namespace
{
// gold vector with nnz for each row
const std::vector<int> goldNNZPerRow =
  {2,3,3,3,2,
   3,4,4,4,3,
   3,4,4,4,3,
   2,3,3,3,2};

// gold vector with correct col ind per row
const std::vector<std::array<std::size_t,4>> goldColIndPerRow =
  {
    {0,1}, {2,1,3}, {4,3,5}, {6,5,7}, {8,7},
    {9, 0, 10}, {11,2,10,12}, {13,4,12,14}, {15,6,14,16}, {17,8,16},
    {18, 9, 19}, {20, 11, 19,21}, {22, 13, 21,23}, {24, 15, 23, 25}, {26, 17, 25},
    {18, 27}, {20, 27,28}, {22, 28,29}, {24,29,30}, {26,30}
  };

// gold vector with correct vals per row
const std::vector< std::array<double,4> > goldValPerRow =
  {
    {0., 0.,},
    {2.0754064337599448e-06, -7.8517110556081196e-08, 6.5322975423424218e-07},
    {2.0754064337599448e-06, -3.6587343239516172e-07, 3.6587343239516172e-07},
    {2.0754064337599448e-06, -6.5322975423424218e-07, 7.8517110556081302e-08},
    {0.,0.},
    {0.,0.},
    {1.3752622897177191e-06, -7.0014414404222578e-07, -6.1489658255606281e-08, 5.1156842203926846e-07},
    {1.3752622897177191e-06, -7.0014414404222578e-07, -2.8652904014743738e-07, 2.8652904014743738e-07},
    {1.3752622897177191e-06, -7.0014414404222578e-07, -5.1156842203926836e-07, 6.148965825560636e-08},
    {0.,0.},
    {0.,0.},
    {1.3151042771684694e-06, -7.6030215659147531e-07, -5.0531292948186884e-08, 4.2039937333284984e-07},
    {1.3151042771684694e-06, -7.6030215659147531e-07, -2.3546533314051836e-07, 2.3546533314051836e-07},
    {1.3151042771684694e-06, -7.6030215659147531e-07, -4.2039937333284973e-07, 5.0531292948186963e-08},
    {0.,0.},
    {0.,0.},
    {-2.0754064337599448e-06, -4.2888015183670156e-08, 3.5681047633576567e-07},
    {-2.0754064337599448e-06, -1.998492457597179e-07, 1.998492457597179e-07},
    {-2.0754064337599448e-06, -3.5681047633576556e-07, 4.2888015183670209e-08},
    {0.,0.}
  };
}

int main(int argc, char *argv[])
{
  Kokkos::initialize (argc, argv);
  {
    std::string sentinel = "PASS";

    using sc_t = double;
    using parser_t = kokkosapp::commonTypes::parser_t;
    using prob_t = kokkosapp::FomProblemRankOneForcing;
    using app_t  = typename prob_t::fom_t;
    using mesh_info_t = typename prob_t::mesh_info_t;
    using forcing_t = typename prob_t::forcing_t;

    parser_t parser(argc, argv);
    mesh_info_t meshInfo(parser.getMeshDir());
    app_t appObj(meshInfo);

    auto matObj = createMaterialModel<sc_t>(parser);
    appObj.computeJacobians(*matObj);

    auto Jd = appObj.viewJacobianDevice(dofId::vp);

    // *** do check ***
    if (Jd.numRows() != 20){ sentinel = "FAIL"; }
    if (Jd.numCols() != 31){ sentinel = "FAIL"; }

    for (std::size_t i=0; i<Jd.numRows(); ++i)
    {
      auto rowview = Jd.rowConst(i);
      const auto l = rowview.length;
      std::cout << "row = " << i << " l = " << l << std::endl;

      // check correct number of nnz for this row
      if (l != goldNNZPerRow[i]){ sentinel = "FAIL"; break; }

      // check that current row has correct col indices
      const int thisRowNNZ = goldNNZPerRow[i];
      const auto goldColInd = goldColIndPerRow[i];
      for (int j=0; j<thisRowNNZ; ++j){
	if (goldColInd[j] != rowview.colidx(j)){ sentinel = "FAIL"; break; }
      }

      // check that current row has correct values
      const auto goldVals = goldValPerRow[i];
      for (int j=0; j<thisRowNNZ; ++j)
      {
	const auto err = std::abs(goldVals[j] - rowview.value(j));
	std::cout << "gold  = " << goldVals[j] << " "
		  << "found = " << rowview.value(j) << " "
		  << "diff = " << err << std::endl;

	if (err > 1e-13){ sentinel = "FAIL"; }
      }
    }

    std::puts(sentinel.c_str());
  }
  Kokkos::finalize();

  return 0;
}
