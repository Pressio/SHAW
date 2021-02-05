
#include "../shared/constants.hpp"
#include "../shared/meta/meta_eigen.hpp"
#include "../shared/meta/meta_kokkos.hpp"
#include "../shared/io/matrix_write.hpp"
#include "../shared/io/matrix_read.hpp"
#include "../shared/io/read_basis.hpp"
#include "../shared/io/read_reference_state.hpp"
#include "../shared/io/vector_write.hpp"
#include "../shared/io/vector_read.hpp"

namespace
{

constexpr auto edyn	= Eigen::Dynamic;
using sc_t		= double;
using snap_t	        = Eigen::Matrix<sc_t, edyn, edyn, Eigen::ColMajor>;
using int_t		= typename snap_t::Index;
using int_p_t		= std::pair<int_t, int_t>;

void loadTargetSnapshotMatrix(const std::string file,
			      int_t & numRows,
			      int_t & numCols,
			      snap_t & M,
			      int useBinary)
{
  // load snapshot matrix
  if (useBinary == 1){
    std::cout << "Using Binary " << std::endl;
    readBinaryMatrixWithSize(file, M);
  }
  else{
    std::cout << "Using ascii " << std::endl;
    readAsciiMatrixWithSize(file, M);
  }
  std::cout << "Done reading file" << std::endl;
  numRows = M.rows();
  numCols = M.cols();
  std::cout << "Current snapshots size: " << numRows << " " << numCols << std::endl;
}

}//end anonym namespace


int main(int argc, char *argv[])
{
  std::string file1 = {};
  std::string file2 = {};
  if(argc<3){
    std::cout << "Not enough args.\n" << std::endl;
    std::cout << "args: file1 file2 tolerance(optional) \n" << std::endl;
    return 1;
  }
  file1 = std::string(argv[1]);
  file2 = std::string(argv[2]);
  std::cout << "Files are: " << file1 << " " << file2 << std::endl;

  // default value of tolerance is 1e-14
  const double defaultTol = 1e-14;
  double tol = defaultTol;
  if (argc==4) tol = std::atof(argv[3]);
  std::cout << "Tolerance = " << tol << std::endl;

  // load snapshot matrix A1
  int_t numRows1  = {}; int_t numSnaps1 = {}; snap_t A1;
  loadTargetSnapshotMatrix(file1, numRows1, numSnaps1, A1, 0);

  // load snapshot matrix A2
  int_t numRows2  = {}; int_t numSnaps2 = {}; snap_t A2;
  loadTargetSnapshotMatrix(file2, numRows2, numSnaps2, A2, 0);

  if (numRows1 != numRows2)
    throw std::runtime_error("numRows1 != numRows2");
  if (numSnaps1 != numSnaps2)
    throw std::runtime_error("numSnaps1 != numSnaps2");

  for (auto j=0; j<numSnaps1; ++j)
  {
    for (auto i=0; i<numRows1; ++i)
    {
      if (std::isnan(A1(i,j)) or std::isnan(A2(i,j)))
	throw std::runtime_error("Found NaN");

      const auto diff = A1(i,j) - A2(i,j);
      if ( std::abs(diff) > tol ){
	std::cout << i << " "
		  << j << " "
		  << A1(i,j) << " "
		  << A2(i,j) << std::endl;
  	throw std::runtime_error("A1(i,j) != A2(i,j)");
      }
    }
  }

  std::cout << "Snaps match!" << std::endl;

  return 0;
}
