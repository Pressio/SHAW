
#include "./fom/fom_problem_rank_one_forcing.hpp"

namespace
{
// gold vector with correct col ind per row
const std::vector<std::array<std::size_t,2>> goldColIndPerRow =
  {
    {5, 0}, {0, 1},{6, 1},{1, 2},{7, 2},{2, 3},{8, 3},{3, 4},{9, 4},
    {10, 5},{5, 6},{11, 6},{6, 7},{12, 7},{7, 8},{13, 8},{8, 9},{14, 9},
    {15, 10},{10, 11},{16, 11},{11, 12},{17, 12},{12, 13},{18, 13},{13, 14},
    {19, 14},{15, 16},{16, 17},{17, 18},{18, 19}
  };

// gold vector with correct vals per row
const std::vector< std::array<double,2> > goldValPerRow =
  {
    {9.1149901848697249e-07,  -1.1639074152729723e-06},
    {-7.127431971039397e-07,  1.9003667686383708e-08},
    {9.1149901848697249e-07,  -1.1639074152729723e-06},
    {-4.2538687526485933e-07,  3.063599895254641e-07},
    {9.1149901848697249e-07,  -1.1639074152729723e-06},
    {-3.0635998952546415e-07,  4.2538687526485928e-07},
    {9.1149901848697249e-07,  -1.1639074152729723e-06},
    {-1.9003667686383771e-08,  7.127431971039397e-07},
    {9.1149901848697249e-07,  -1.1639074152729723e-06},
    {9.361906800816778e-07,  -1.139215753678267e-06},
    {-5.5817560406309572e-07,  1.4882476231779004e-08},
    {9.361906800816778e-07,  -1.139215753678267e-06},
    {-3.3313622217126481e-07,  2.3992185812361002e-07},
    {9.361906800816778e-07,  -1.139215753678267e-06},
    {-2.3992185812361007e-07,  3.3313622217126475e-07},
    {9.361906800816778e-07,  -1.139215753678267e-06},
    {-1.4882476231779053e-08,  5.5817560406309572e-07},
    {9.361906800816778e-07,  -1.139215753678267e-06},
    {9.5280156129768852e-07,  -1.1226048724622563e-06},
    {-4.5870046712890708e-07,  1.2230199152129572e-08},
    {9.5280156129768852e-07,  -1.1226048724622563e-06},
    {-2.7376642693657577e-07,  1.9716423934446095e-07},
    {9.5280156129768852e-07,  -1.1226048724622563e-06},
    {-1.9716423934446101e-07,  2.7376642693657571e-07},
    {9.5280156129768852e-07,  -1.1226048724622563e-06},
    {-1.2230199152129614e-08,  4.5870046712890708e-07},
    {9.5280156129768852e-07,  -1.1226048724622563e-06},
    {-3.8931821157145032e-07,  1.038027994798545e-08},
    {-2.3235698099540268e-07,  1.6734151052403315e-07},
    {-1.6734151052403317e-07,  2.3235698099540265e-07},
    {-1.0380279947985484e-08,  3.8931821157145032e-07}
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

    auto matObj = createMaterialModel<sc_t>(parser, meshInfo);
    appObj.computeJacobians(*matObj);

    auto Jd = appObj.viewJacobianDevice(dofId::sp);

    // *** do check ***
    if (Jd.numRows() != 31){ sentinel = "FAIL"; }
    if (Jd.numCols() != 20){ sentinel = "FAIL"; }

    for (std::size_t i=0; i<Jd.numRows(); ++i)
    {
      auto rowview = Jd.rowConst(i);
      const auto l = rowview.length;
      std::cout << "row = " << i << " l = " << l << std::endl;

      // check correct number of nnz for this row
      if (l != 2){ sentinel = "FAIL"; break; }

      // check that current row has correct col indices
      const auto goldColInd = goldColIndPerRow[i];
      //std::cout << rowview.colidx(0) << " " << rowview.colidx(1) << std::endl;
      if (goldColInd[0] != rowview.colidx(0)){ sentinel = "FAIL"; break; }
      if (goldColInd[1] != rowview.colidx(1)){ sentinel = "FAIL"; break; }

      // check that current row has correct values
      const auto goldVals = goldValPerRow[i];
      for (int j=0; j<2; ++j)
      {
	const auto err = std::abs(goldVals[j] - rowview.value(j));
	std::cout << "gold = " << goldVals[j] <<  " "
		  << "found = " << rowview.value(j) <<  " "
		  << "diff = " << err << std::endl;

	if (err > 1e-13){ sentinel = "FAIL"; }
      }
    }

    std::puts(sentinel.c_str());
  }
  Kokkos::finalize();

  return 0;
}
