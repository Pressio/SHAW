
#include "./shared/all.hpp"
#include "./kokkos/common_types.hpp"
#include "./kokkos/shwavepp.hpp"

template <typename T>
typename std::enable_if<
  !std::is_same<typename T::memory_space, Kokkos::HostSpace>::value
  >::type
runTest(T J, std::string & sentinel)
{
  const std::vector<std::size_t> goldInd =
    {
     0,1, 2,1,3, 4,3,5, 6,5,7, 8,7,
     9, 0, 10, 11,2,10,12, 13,4,12,14, 15,6,14,16, 17,8,16,
     18, 9, 19, 20, 11, 19,21, 22, 13, 21,23, 24, 15, 23, 25, 26, 17, 25,
     18, 27, 20, 27,28, 22, 28,29, 24,29,30, 26,30
    };

  const std::vector<int> goldPtr =
    {
     0,
     2,
     5,
     8,
     11,
     13,
     16,
     20,
     24,
     28,
     31,
     34,
     38,
     42,
     46,
     49,
     51,
     54,
     57,
     60,
     62
    };

  const std::vector<double> goldValues =
    {
     0., 0.,
     2.0754064337599448e-06, -7.8517110556081196e-08, 6.5322975423424218e-07,
     2.0754064337599448e-06, -3.6587343239516172e-07, 3.6587343239516172e-07,
     2.0754064337599448e-06, -6.5322975423424218e-07, 7.8517110556081302e-08,
     0.,0.,
     0.,0.,0.,
     1.3752622897177191e-06, -7.0014414404222578e-07, -6.1489658255606281e-08, 5.1156842203926846e-07,
     1.3752622897177191e-06, -7.0014414404222578e-07, -2.8652904014743738e-07, 2.8652904014743738e-07,
     1.3752622897177191e-06, -7.0014414404222578e-07, -5.1156842203926836e-07, 6.148965825560636e-08,
     0.,0.,0.,
     0.,0.,0.,
     1.3151042771684694e-06, -7.6030215659147531e-07, -5.0531292948186884e-08, 4.2039937333284984e-07,
     1.3151042771684694e-06, -7.6030215659147531e-07, -2.3546533314051836e-07, 2.3546533314051836e-07,
     1.3151042771684694e-06, -7.6030215659147531e-07, -4.2039937333284973e-07, 5.0531292948186963e-08,
     0.,0.,0.,
     0.,0.,
     -2.0754064337599448e-06, -4.2888015183670156e-08, 3.5681047633576567e-07,
     -2.0754064337599448e-06, -1.998492457597179e-07, 1.998492457597179e-07,
     -2.0754064337599448e-06, -3.5681047633576556e-07, 4.2888015183670209e-08,
     0.,0.
    };

  using jd_t = decltype(J);

  using rm_d_t = typename jd_t::row_map_type::non_const_type;
  using ind_d_t = typename jd_t::index_type::non_const_type;
  using v_d_t = typename jd_t::values_type;
  auto val_h = Kokkos::create_mirror_view (J.values);
  auto ptr_h = Kokkos::create_mirror_view (J.graph.row_map);
  auto ind_h = Kokkos::create_mirror_view (J.graph.entries);
  Kokkos::deep_copy(val_h, J.values);
  Kokkos::deep_copy(ptr_h, J.graph.row_map);
  Kokkos::deep_copy(ind_h, J.graph.entries);

  for (std::size_t i=0; i<val_h.extent(0); ++i)
  {
    const auto err1 = std::abs(goldValues[i] - val_h(i));
    const auto err2 = (double) goldInd[i] - (double) ind_h(i);
    if (err1 > 1e-13){ sentinel = "FAIL"; }
    if (err2 > 1e-13){ sentinel = "FAIL"; }
    //std::cout << ind_h(i) << " " << val_h(i) << "\n";
  }

  for (std::size_t i=0; i<ptr_h.extent(0); ++i)
  {
    const auto err2 = (double) goldPtr[i] - (double) ptr_h(i);
    if (err2 > 1e-13){ sentinel = "FAIL"; }
  }
}

template <typename T>
typename std::enable_if<
  std::is_same<typename T::memory_space, Kokkos::HostSpace>::value
  >::type
runTest(T J, std::string & sentinel)
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
     {0.,0.,0.},
     {1.3752622897177191e-06, -7.0014414404222578e-07, -6.1489658255606281e-08, 5.1156842203926846e-07},
     {1.3752622897177191e-06, -7.0014414404222578e-07, -2.8652904014743738e-07, 2.8652904014743738e-07},
     {1.3752622897177191e-06, -7.0014414404222578e-07, -5.1156842203926836e-07, 6.148965825560636e-08},
     {0.,0.,0.},
     {0.,0.,0.},
     {1.3151042771684694e-06, -7.6030215659147531e-07, -5.0531292948186884e-08, 4.2039937333284984e-07},
     {1.3151042771684694e-06, -7.6030215659147531e-07, -2.3546533314051836e-07, 2.3546533314051836e-07},
     {1.3151042771684694e-06, -7.6030215659147531e-07, -4.2039937333284973e-07, 5.0531292948186963e-08},
     {0.,0.,0.},
     {0.,0.},
     {-2.0754064337599448e-06, -4.2888015183670156e-08, 3.5681047633576567e-07},
     {-2.0754064337599448e-06, -1.998492457597179e-07, 1.998492457597179e-07},
     {-2.0754064337599448e-06, -3.5681047633576556e-07, 4.2888015183670209e-08},
     {0.,0.}
    };

  // *** do check ***
  if (J.numRows() != 20){ sentinel = "FAIL"; }
  if (J.numCols() != 31){ sentinel = "FAIL"; }

  for (std::size_t i=0; i<J.numRows(); ++i)
    {
      auto rowview = J.rowConst(i);

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
}

int main(int argc, char *argv[])
{
  Kokkos::initialize (argc, argv);
  {
    std::string sentinel = "PASS";

    using sc_t = typename kokkosapp::commonTypes::scalar_type;
    using parser_t = typename kokkosapp::commonTypes::parser_type;
    using mesh_info_t = typename kokkosapp::commonTypes::mesh_info_type;

    parser_t parser(argc, argv);
    mesh_info_t meshInfo(parser.getMeshDir());

    auto matObj = createMaterialModel<sc_t>(parser, meshInfo);
    kokkosapp::ShWavePP<kokkosapp::commonTypes> appObj(meshInfo);
    appObj.computeJacobians(*matObj);

    auto Jd = appObj.viewJacobianDevice(dofId::vp);
    runTest(Jd, sentinel);
    std::puts(sentinel.c_str());
  }
  Kokkos::finalize();

  return 0;
}
