
#include "CLI11.hpp"
#include "Eigen/Dense"
#include "../shared/constants.hpp"

namespace
{

template <typename mat_t>
void write_matrix_to_ascii(const std::string fileName,
			   const mat_t & A,
			   std::size_t m,
			   std::size_t n,
			   bool writeSize = true)
{
  std::ofstream file; file.open(fileName);
  if (writeSize){
    file << m << " " << n << std::endl;
  }

  for (std::size_t i=0; i<m; i++){
    for (std::size_t j=0; j<n; j++){
      file << std::setprecision(15) << A(i,j) << " ";
    }
    file << std::endl;
  }
  file.close();
}

template<typename sc_t>
void write_contig_matrix_to_binary(const std::string filename,
				   const sc_t * A,
				   std::size_t m,
				   std::size_t n,
				   bool writeSize = true)
{
  std::ofstream out(filename, std::ios::out | std::ios::binary | std::ios::trunc);
  if (writeSize){
    out.write((char*) (&m), sizeof(std::size_t));
    out.write((char*) (&n), sizeof(std::size_t));
  }
  out.write((char*) A, m*n*sizeof(sc_t) );
  out.close();
}

using sc_t     = double;
using snap_t   = Eigen::Matrix<sc_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
using vec_t    = Eigen::Matrix<sc_t, Eigen::Dynamic, 1>;
using int_t    = typename snap_t::Index;
using int_p_t  = std::pair<int_t, int_t>;

template<class T>
void writeToFile(const std::string fileName,
		 const T & A,
		 const std::string & outformat,
		 const bool writeSize = true)
{
  write_contig_matrix_to_binary(fileName, A.data(),
				A.rows(), A.cols(), writeSize);
}

template<class dmat_t>
void readBinaryMatrixWithSize(const std::string filename, dmat_t & M)
{
  using int_t = typename dmat_t::Index;
  using sc_t  = typename dmat_t::Scalar;
  std::ifstream fin(filename, std::ios::in | std::ios::binary);
  fin.exceptions(std::ifstream::failbit | std::ifstream::badbit);

  std::size_t rows={};
  std::size_t cols={};
  fin.read((char*) (&rows),sizeof(std::size_t));
  fin.read((char*) (&cols),sizeof(std::size_t));
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

snap_t loadTargetSnapshotMatrix(const std::string file,
			      std::size_t & numRows,
			      std::size_t & numCols,
			      const std::string & inputFormat)
{
  snap_t M;

  if (inputFormat == "binary"){
    std::cout << "Using Binary " << std::endl;
    readBinaryMatrixWithSize(file, M);
  }
  // else{
  //   std::cout << "Using ascii " << std::endl;
  //   fillMatrixFromAscii(file, M);
  // }

  numRows = M.rows();
  numCols = M.cols();
  std::cout << "Current snapshots size: "
	    << numRows << " "
	    << numCols << std::endl;
  return M;
}

template <typename T>
void doThinSVDm0(T & A,
		 const std::string lsvFileName,
		 const std::string rsvFileName,
		 const std::string singValuesFileName,
		 const std::string dofType,
		 const std::string outputFormat)
{

  std::cout << "Computing SVD method 1" << std::endl;

  std::cout << "MinMax(A) "
	    << A.minCoeff() << " "
	    << A.maxCoeff() << std::endl;

  //Eigen::JacobiSVD<T> svd(A, Eigen::ComputeThinU);
  Eigen::BDCSVD<T> svd(A, Eigen::ComputeThinU);
  std::cout << "svd_matrix_" + dofType + "_rank = " << svd.rank() << std::endl;

  // get sing values and left sing vectors
  const auto & singVal = svd.singularValues();
  const auto & U       = svd.matrixU();

  {
    std::cout << "Printing sing values" << std::endl;
    std::ofstream file; file.open(singValuesFileName);
    file << std::setprecision(dblFmt) << singVal << '\n';
    file.close();
  }

  std::cout << std::endl;
  std::cout << "Printing left-sing vectors to file" << std::endl;
  {
    writeToFile(lsvFileName, U, outputFormat);
  }

  std::cout << "Done with SVD" << std::endl;
  std::cout << "----------------" << std::endl;
}


template <typename T>
void doThinSVDm1(T & A,
		 const std::string lsvFileName,
		 const std::string rsvFileName,
		 const std::string singValuesFileName,
		 const std::string dofType,
		 const std::string outputFormat)
{

  std::cout << "Computing SVD method 2" << std::endl;

  std::cout << "MinMax(A) "
	    << A.minCoeff() << " "
	    << A.maxCoeff() << std::endl;

  Eigen::HouseholderQR<snap_t> qr(A);

  using R_type = Eigen::Matrix<sc_t, -1,-1>;
  R_type R = qr.matrixQR().block(0,0,A.cols(),A.cols()).template triangularView<Eigen::Upper>();

  Eigen::BDCSVD<R_type> svd(R, Eigen::ComputeFullU);
  const auto singVal = svd.singularValues();
  const auto U0      = svd.matrixU();

  snap_t thinQ = snap_t::Identity(A.rows(),A.cols());
  thinQ = qr.householderQ() * thinQ;
  snap_t U = thinQ * U0;

  {
    std::cout << "Printing sing values" << std::endl;
    std::ofstream file; file.open(singValuesFileName);
    file << std::setprecision(dblFmt) << singVal << '\n';
    file.close();
  }

  std::cout << std::endl;
  std::cout << "Printing left-sing vectors to file" << std::endl;
  {
    writeToFile(lsvFileName, U, outputFormat);
  }

  std::cout << "Done with SVD" << std::endl;
  std::cout << "----------------" << std::endl;
}

void processDirs(std::string dofName,
		 const std::vector<std::string> & dirs,
		 const std::string & outputFormat,
		 const std::string & inputFormat,
		 const int method)
{
  // the matrix containing all snapshots for this dof
  snap_t allSnaps;

  // loop over dirs
  auto startTime1 = std::chrono::high_resolution_clock::now();
  for (std::size_t iDir=0; iDir<dirs.size(); ++iDir)
  {
    const auto dirName = dirs[iDir];
    const std::string file = dirName + "/snaps_"+dofName+"_0";
    std::cout << "\nReading snapshots: " << file << std::endl;

    // load curren snapshot matrix
    std::size_t numRows  = {};
    std::size_t numSnaps = {};
    auto snaps = loadTargetSnapshotMatrix(file, numRows, numSnaps, inputFormat);

    if (iDir == 0){
      // when iDir == 0, needs to full resize allSnaps
      allSnaps.resize(numRows, numSnaps);
      allSnaps = snaps;
    }
    else{
      // we need to increase size to make space for new snapshots but
      // making sure the currently read snapshot matrix has
      // the same number of rows otherwise something is wrong.

      // current number of rows in allSnaps
      const auto currNRowsAllSnaps = allSnaps.rows();
      if (numRows != currNRowsAllSnaps){
	throw std::runtime_error("Mismatching # rows of current data with previous");
      }

      // current number of cols in allSnaps
      const auto currNColsAllSnaps = allSnaps.cols();
      // resize to make space for current data
      allSnaps.conservativeResize(numRows, currNColsAllSnaps+numSnaps);
      // copy current snapshot matrix into allSnaps
      allSnaps.block(0, currNColsAllSnaps, numRows, numSnaps) = snaps;
    }
  }//loop over dirs

  std::cout << "Final snapshot matrix size: "
	    << allSnaps.rows() << " "
	    << allSnaps.cols()  << std::endl;

  const auto finishTime1 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed1 = finishTime1 - startTime1;
  std::cout << "svd:read:" << dofName << ": "
	    << std::fixed << std::setprecision(10)
	    << elapsed1.count() << std::endl;

  auto startTime2 = std::chrono::high_resolution_clock::now();
  auto lsvFN = "lsv_"+dofName;
  auto rsvFN = "rsv_"+dofName;
  auto svFN  = "sva_"+dofName;

  if (method==1){
    doThinSVDm0(allSnaps, lsvFN, rsvFN, svFN, dofName, outputFormat);
  }
  else if (method==2){
    doThinSVDm1(allSnaps, lsvFN, rsvFN, svFN, dofName, outputFormat);
  }

  const auto finishTime2 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed2 = finishTime2 - startTime1;
  std::cout << "svd:time" << dofName << ": "
	    << std::fixed << std::setprecision(10)
	    << elapsed2.count() << std::endl;

}//end processDirs

}//end anonym namespace

int main(int argc, char *argv[])
{
  CLI::App app;
  std::vector<std::string> dirs = {};
  std::string outputFormat = {};
  std::string inputFormat = {};
  int method = -1;
  app.add_option("--dirs", dirs,
		 "Directories")->required();
  app.add_option("--informat", inputFormat,
		 "Inputformat: binary/ascii")->required();
  app.add_option("--outformat", outputFormat,
		 "Outputformat: binary/ascii")->required();
  app.add_option("--method", method,
		 "Method: 1: Eigen SVD, 2: via QR")->required();
  try{
    app.parse(argc, argv);
  }
  catch (...){}

  if (method == -1){
    throw std::runtime_error("Invalid method");
  }

  processDirs("vp", dirs, inputFormat, outputFormat, method);
  processDirs("sp", dirs, inputFormat, outputFormat, method);

  return 0;
}



// template <typename T>
// void doThinSVDm1(T & S,
// 		 const std::string lsvFileName,
// 		 const std::string rsvFileName,
// 		 const std::string singValuesFileName,
// 		 const std::string dofType,
// 		 const int useBinary,
// 		 const int doScaling)
// {

//   std::cout << "Computing SVD method 1" << std::endl;

//   using mat_t = Eigen::Matrix<sc_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;

//   std::cout << "MinMax(A) "
// 	    << S.minCoeff() << " "
// 	    << S.maxCoeff() << std::endl;

//   // compute B = S^T S
//   mat_t B = S.transpose() * S;
//   const auto maxValue = B.maxCoeff();
//   const auto minValue = B.minCoeff();
//   //if (doScaling==1) scaleSnapshots(B, minValue, maxValue);

//   std::cout << "MinMax(B) "
// 	    << B.minCoeff() << " "
// 	    << B.maxCoeff() << std::endl;

//   // do svd on B
//   Eigen::BDCSVD<T> svd(B, Eigen::ComputeThinU);
//   svd.setThreshold(1e-16);
//   const auto rankk = svd.rank();
//   std::cout << "svd_matrix_" + dofType + "_rank = " << rankk << std::endl;

//   // get sing values and left sing vectors
//   const auto & singVal = svd.singularValues();
//   const auto & V       = svd.matrixU();
//   //if (doScaling==1) unscale(V, minValue, maxValue);

//   Eigen::VectorXd sv2(singVal.size());
//   Eigen::VectorXd sv3(singVal.size());
//   for (auto i=0; i<sv2.size(); ++i){
//     sv2(i) = 1./std::sqrt(singVal(i));
//     sv3(i) = std::sqrt(singVal(i));
//   }

//   {
//     std::cout << "Printing sing values" << std::endl;
//     std::ofstream file; file.open(singValuesFileName);
//     file << std::setprecision(dblFmt) << sv3 << '\n';
//     file.close();
//   }

//   mat_t U = S * V * sv2.asDiagonal();

//   std::cout << std::endl;
//   std::cout << "Printing left-sing vectors to file" << std::endl;
//   {
//     writeToFile(lsvFileName, U, useBinary);
//   }

//   std::cout << "Done with SVD" << std::endl;
//   std::cout << "----------------" << std::endl;
// }





  // using T = Eigen::MatrixXd;
  // T A = T::Random(20,4);
  // write_matrix_to_ascii("A.txt", A, A.rows(), A.cols(), false);

  // {
  //   Eigen::BDCSVD<T> svd(A, Eigen::ComputeThinU);
  //   const auto & singVal = svd.singularValues();
  //   const auto & U       = svd.matrixU();
  //   std::cout << singVal << '\n';
  //   std::cout << "\n";
  //   std::cout << U << '\n';
  // }

  // {
  //   int m = A.cols();
  //   Eigen::HouseholderQR<T> qr(A);
  //   T Q = qr.householderQ();
  //   T R = qr.matrixQR().block(0,0,m,m).template triangularView<Eigen::Upper>();
  //   Eigen::BDCSVD<decltype(R)> svd(R, Eigen::ComputeFullU | Eigen::ComputeFullV);
  //   const auto & singVal = svd.singularValues();
  //   const auto & U0      = svd.matrixU();
  //   const auto & V0      = svd.matrixV();

  //   T U = Q * U0;
  //   std::cout << "\n";
  //   std::cout << singVal << '\n';
  //   std::cout << "\n";
  //   std::cout << U << '\n';
  //   std::cout << U.col(0).dot(U.col(1)) << '\n';
  // }