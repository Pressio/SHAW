
#include <fstream>
#include <iomanip>
#include <iostream>
#include <ostream>
#include "../shared/constants.hpp"
#include "../shared/meta/meta_eigen.hpp"
// #include "../shared/io/matrix_read.hpp"
// #include "../shared/io/vector_read.hpp"

namespace
{

template<class dmat_t>
typename std::enable_if< is_dynamic_matrix_eigen<dmat_t>::value >::type
readBinaryMatrixWithSize(const std::string filename, dmat_t & M)
{
  using int_t = typename dmat_t::Index;
  using sc_t  = typename dmat_t::Scalar;
  std::ifstream fin(filename, std::ios::in | std::ios::binary);
  fin.exceptions(std::ifstream::failbit | std::ifstream::badbit);

  int_t rows={};
  int_t cols={};
  fin.read((char*) (&rows),sizeof(int_t));
  fin.read((char*) (&cols),sizeof(int_t));
  const auto nBytes = rows*cols*sizeof(sc_t);
  M.resize(rows, cols);
  fin.read( (char *) M.data(), nBytes );

  if (!fin){
    std::cout << std::strerror(errno) << std::endl;
    throw std::runtime_error("ERROR READING binary file");
  }
  else
    std::cout << fin.gcount() << " bytes read\n";

  fin.close();
}

template <typename dmat_t>
typename std::enable_if<is_dynamic_matrix_eigen<dmat_t>::value>::type
readAsciiMatrixWithSize(const std::string fileName, dmat_t & M)
{
  using int_t = typename dmat_t::Index;
  std::ifstream source; source.open(fileName, std::ios_base::in);
  std::string line, colv;

  {
    // first line contains the size of the matrix
    std::getline(source, line);
    std::istringstream in(line);
    std::string col1, col2;
    in >> col1; in >> col2;
    M.resize( std::stoi(col1), std::stoi(col2) );
  }

  // then read the actual data
  int_t iRow = 0;
  while (std::getline(source, line) ){
    std::istringstream in(line);
    for (int_t j=0; j<M.cols(); ++j){
      in >> colv;
      M(iRow, j) = atof(colv.c_str());
    }
    iRow++;
  }
  source.close();
}

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
  	throw std::runtime_error("Wrong matching A1(i,j) != A2(i,j)");
      }
    }
  }

  std::cout << "Snaps match!" << std::endl;

  return 0;
}
