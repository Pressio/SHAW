
#include "../shared/io.hpp"
#include "../shared/constants.hpp"
#include "../shared/metaf.hpp"

// namespace
// {
// using sc_t		= double;

// constexpr auto edy	= Eigen::Dynamic;
// using snap_t	        = Eigen::Matrix<sc_t, edy, edy, Eigen::ColMajor>;
// using vec_t	        = Eigen::Matrix<sc_t, edy, 1>;
// using int_t		= typename snap_t::StorageIndex;

// using int_p_t		= std::pair<int_t, int_t>;
// using dir_list_t	= std::vector<std::string>;

//static_assert( std::is_same< typename snap_t::StorageIndex, int32_t >::value, "" );

// void loadTargetSnapshotMatrix(const std::string file,
// 			      int_t & numRows,
// 			      int_t & numCols,
// 			      snap_t & M,
// 			      int useBinary)
// {
//   // load snapshot matrix
//   if (useBinary == 1){
//     std::cout << "Using Binary " << std::endl;
//     readBinaryMatrixWithSize(file, M);
//   }
//   else{
//     std::cout << "Using ascii " << std::endl;
//     readAsciiMatrixWithSize(file, M);
//   }
//   numRows = M.rows();
//   numCols = M.cols();
//   std::cout << "Current snapshots size: " << numRows << " " << numCols << std::endl;
// }
// }

// void doMatrixTest(int_t nR, int_t nC, int useBinary)
// {
//   const std::string fileName = "matFile";
//   snap_t M(nR, nC);
//   //M.setRandom();
//   snap_t M1;

//   std::cout << "Printing matrix to file" << std::endl;
//   writeMatrixWithSizeToFile(fileName, M, useBinary);
//   std::cout << "Done printing matrix to file\n" << std::endl;

//   std::cout << "Reading matrix from file" << std::endl;
//   loadTargetSnapshotMatrix(fileName, nR, nC, M1, useBinary);
//   std::cout << "Done reading matrix to file" << std::endl;

//   if (M != M1){
//     throw std::runtime_error("Matrices don't match");
//   }
// }


template<class T>
void writeStdVectorBinary(const std::string filename, const T& a)
{
  using sc_t  = typename T::value_type;
  std::ofstream out(filename, std::ios::out | std::ios::binary | std::ios::trunc);
  std::size_t rows=a.size();
  //out.write((char*) (&rows), sizeof(std::size_t));
  out.write((char*) a.data(), rows*sizeof(sc_t) );
  out.close();
}

template<class vec_t>
void readBinaryStdVector(std::string filename, vec_t & a, std::size_t fs)
{
  using sc_t  = typename vec_t::value_type;

  std::ifstream in(filename, std::ios::in | std::ios::binary);
  in.exceptions(std::ifstream::failbit | std::ifstream::badbit);

  std::size_t rows= fs/sizeof(sc_t);
  // in.read((char*) (&rows),sizeof(std::size_t));
  std::size_t nBytes = rows*sizeof(sc_t);
  std::cout << "Bytes to read " << nBytes << "\n";

  in.read( (char *) &(*a.data()), nBytes );

  if (!in)
    std::cout << "ERROR READING" << std::strerror(errno) << std::endl;
  else
    std::cout << in.gcount() << " bytes read\n";

  in.close();
  std::cout << std::setprecision(18) << a[0] << " " << a[rows-1] << std::endl;
}

std::size_t printSize(const std::string& address) {
  std::fstream motd(address.c_str(), std::ios::binary|std::ios::in|std::ios::ate);
  if(motd) {
    std::fstream::pos_type size = motd.tellg();
    std::cout << address << " " << size << "\n";
    return size;
  } else {
    perror(address.c_str());
  }
}

void doVectorTest(std::size_t nR, std::size_t nC)
{
  const std::string fileName = "matFile";

  // create vector
  std::size_t totN = nR*nC;
  Eigen::VectorXd vec(totN);
  vec.setRandom();
  std::cout << std::setprecision(18)
	    << vec[0] << " " << vec[totN-1] << std::endl;

  std::cout << "Printing vector to file" << std::endl;
  writeStdVectorBinary(fileName, vec);
  std::cout << "Done printing vector to file\n" << std::endl;

  std::cout << "Calculating size" << std::endl;
  const std::size_t fsize = printSize(fileName);

  std::cout << "Reading vector from file" << std::endl;
  Eigen::VectorXd vec2(fsize);
  readBinaryStdVector(fileName, vec2, fsize);
  std::cout << "Done reading vector to file" << std::endl;
}

int main(int argc, char *argv[])
{
  std::size_t nR = {};
  std::size_t nC = {};
  int useBinary   = 0;

  if (argc!=4){
    std::cout << "Not enough args.\n" << std::endl;
    std::cout << "args: nRows nCols useBinary \n" << std::endl;
    return 1;
  }

  nR = std::atoi(argv[1]);
  nC = std::atoi(argv[2]);
  useBinary = std::atoi(argv[3]);
  std::cout << nR << " " << nC << " " << useBinary << std::endl;

  //doMatrixTest(nR, nC, useBinary);
  doVectorTest(nR, nC);

  return 0;
}
