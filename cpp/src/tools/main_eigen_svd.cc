
#include "../shared/constants.hpp"
#include "../shared/meta.hpp"
#include "../shared/io.hpp"

namespace
{

using sc_t		= double;
using snap_t	        = Eigen::Matrix<sc_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
using vec_t	        = Eigen::Matrix<sc_t, Eigen::Dynamic, 1>;
using int_t		= typename snap_t::Index;
using int_p_t		= std::pair<int_t, int_t>;
using dir_list_t	= std::vector<std::string>;

// template <typename T>
// void scaleSnapshots(T & S, sc_t minVal, sc_t maxVal)
// {
//   // constexpr auto one = constants::one<sc_t>();
//   // constexpr auto two = constants::two<sc_t>();
//   // const auto a = two / (maxValue - minValue);
//   // const auto b = one - a * maxValue;
//   for (auto i=0; i<S.rows(); ++i){
//     for (auto j=0; j<S.cols(); j++){
//       S(i,j) = (S(i,j) - minVal)/(maxVal-minVal);
//     }
//   }
//   // for (auto j=0; j<S.cols(); j++)
//   //   std::cout << std::setprecision(dblFmt) << S.col(j).mean() << " ";

//   //std::cout << " " << std::endl;
//   //std::cout << std::setprecision(dblFmt) << a << " " << b << std::endl;
// }

// template <typename T>
// void unscale(T & A, sc_t minVal, sc_t maxVal)
// {
//   for (auto i=0; i<A.rows(); ++i){
//     for (auto j=0; j<A.cols(); j++){
//       A(i,j) = A(i,j)*(maxVal-minVal) + minVal;
//     }
//   }
//   // for (auto j=0; j<S.cols(); j++)
//   //   std::cout << std::setprecision(dblFmt) << S.col(j).mean() << " ";
//   //std::cout << " " << std::endl;
//   //std::cout << std::setprecision(dblFmt) << a << " " << b << std::endl;
// }


// template <typename T>
// void scaleSnapshots2(T & S, const std::string dofType, const int useBinary)
// {
//   // constexpr auto zero = constants::zero<sc_t>();
//   // vec_t refVec(S.rows());
//   // for (auto i=0; i<S.rows(); ++i){
//   //   refVec(i) = zero; //S.row(i).mean();
//   // }

//   // writeVectorWithSizeToFile("ref_state_" + dofType, refVec, useBinary);

//   // for (auto i=0; i<S.rows(); ++i){
//   //   for (auto j=0; j<S.cols(); j++){
//   //     S(i,j) = S(i,j) - refVec(i);
//   //   }
//   // }
//   for (auto j=0; j<S.cols(); j++)
//     std::cout << std::setprecision(dblFmt) << S.col(j).mean() << " ";

//   std::cout << "\n----\n" << std::endl;
//   for (auto j=0; j<S.cols(); j++)
//     std::cout << std::setprecision(dblFmt) << S.col(j).maxCoeff() << " ";

//   std::cout << "\n----\n" << std::endl;
//   for (auto j=0; j<S.cols(); j++)
//     std::cout << std::setprecision(dblFmt) << S.col(j).minCoeff() << " ";

//   // for (auto j=0; j<S.rows(); j++)
//   //   std::cout << std::setprecision(dblFmt) << S.row(j).mean() << " ";
// }

template <typename T>
void doThinSVDm0(T & S,
		 const std::string lsvFileName,
		 const std::string rsvFileName,
		 const std::string singValuesFileName,
		 const std::string dofType,
		 const int useBinary,
		 const int doScaling)
{

  std::cout << "Computing SVD method 0" << std::endl;

  std::cout << "MinMax(A) "
	    << S.minCoeff() << " "
	    << S.maxCoeff() << std::endl;

  //Eigen::JacobiSVD<T> svd(S, Eigen::ComputeThinU);
  Eigen::BDCSVD<T> svd(S, Eigen::ComputeThinU); // | Eigen::ComputeThinV);
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
    writeToFile(lsvFileName, U, useBinary);
  }

  std::cout << "Done with SVD" << std::endl;
  std::cout << "----------------" << std::endl;
}


template <typename T>
void doThinSVDm1(T & S,
		 const std::string lsvFileName,
		 const std::string rsvFileName,
		 const std::string singValuesFileName,
		 const std::string dofType,
		 const int useBinary,
		 const int doScaling)
{

  std::cout << "Computing SVD method 1" << std::endl;

  using mat_t = Eigen::Matrix<sc_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;

  std::cout << "MinMax(A) "
	    << S.minCoeff() << " "
	    << S.maxCoeff() << std::endl;

  // compute B = S^T S
  mat_t B = S.transpose() * S;
  const auto maxValue = B.maxCoeff();
  const auto minValue = B.minCoeff();
  //if (doScaling==1) scaleSnapshots(B, minValue, maxValue);

  std::cout << "MinMax(B) "
	    << B.minCoeff() << " "
	    << B.maxCoeff() << std::endl;

  // do svd on B
  Eigen::BDCSVD<T> svd(B, Eigen::ComputeThinU);
  svd.setThreshold(1e-16);
  const auto rankk = svd.rank();
  std::cout << "svd_matrix_" + dofType + "_rank = " << rankk << std::endl;

  // get sing values and left sing vectors
  const auto & singVal = svd.singularValues();
  const auto & V       = svd.matrixU();
  //if (doScaling==1) unscale(V, minValue, maxValue);

  Eigen::VectorXd sv2(singVal.size());
  Eigen::VectorXd sv3(singVal.size());
  for (auto i=0; i<sv2.size(); ++i){
    sv2(i) = 1./std::sqrt(singVal(i));
    sv3(i) = std::sqrt(singVal(i));
  }

  {
    std::cout << "Printing sing values" << std::endl;
    std::ofstream file; file.open(singValuesFileName);
    file << std::setprecision(dblFmt) << sv3 << '\n';
    file.close();
  }

  mat_t U = S * V * sv2.asDiagonal();

  std::cout << std::endl;
  std::cout << "Printing left-sing vectors to file" << std::endl;
  {
    writeToFile(lsvFileName, U, useBinary);
  }

  std::cout << "Done with SVD" << std::endl;
  std::cout << "----------------" << std::endl;
}


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
  numRows = M.rows();
  numCols = M.cols();
  std::cout << "Current snapshots size: " << numRows << " " << numCols << std::endl;
}


void processDirs(std::string dofName,
		 const dir_list_t & dirs,
		 const int doTranspose,
		 const int useBinary,
		 const int method,
		 const int doScaling)
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
    int_t numRows  = {};
    int_t numSnaps = {};
    snap_t snaps;
    loadTargetSnapshotMatrix(file, numRows, numSnaps, snaps, useBinary);

    // store data
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
      if (numRows != currNRowsAllSnaps)
	throw std::runtime_error("Mismatching # rows of current data with previous");

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
  // file names for:
  // left-sing vectors (lsv), right-sing vectors (rsv) and sing values (sing_vals)
  auto lsvFN = "lsv_"+dofName;
  auto rsvFN = "rsv_"+dofName;
  auto svFN  = "sva_"+dofName;
  // do thin svd
  if (doTranspose == 1){
    std::cout << "Computing SVD on A.transpose()" << std::endl;
    lsvFN += "_T"; rsvFN += "_T"; svFN  += "_T";
    snap_t AT = allSnaps.transpose();
    doThinSVDm0(AT, lsvFN, rsvFN, svFN, dofName, useBinary, doScaling);
  }
  else{
    if(method==0){
      doThinSVDm0(allSnaps, lsvFN, rsvFN, svFN, dofName, useBinary, doScaling);
    }
    else{
      doThinSVDm1(allSnaps, lsvFN, rsvFN, svFN, dofName, useBinary, doScaling);
    }
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
  // *** parse args ***
  // vector of directories to read snapshots from
  dir_list_t dirs;
  auto numDirsIn  = -1;
  int doTranspose = 0;
  int useBinary   = 0;

  // 0: do svd on A
  // 1: use A^T A
  int method = 0;

  int doScaling = 0;

  if(argc>=2){
    numDirsIn = argc-3;
    std::cout << "Directories Passed: " << std::endl;
    for(auto i=0; i<numDirsIn; ++i){
      dirs.push_back( std::string(argv[i+1]) );
      std::cout << dirs[i] << std::endl;
    }

    // use binary or ascii
    useBinary = std::atoi(argv[argc-2]);
    // method
    method = std::atoi(argv[argc-1]);
    // // scaling
    // doScaling = std::atoi(argv[argc-1]);
  }

  processDirs("vp", dirs, doTranspose, useBinary, method, doScaling);
  processDirs("sp", dirs, doTranspose, useBinary, method, doScaling);

  return 0;
}
