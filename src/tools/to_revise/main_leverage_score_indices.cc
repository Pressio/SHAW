
#include "UTILS_ALL"
#include "CONTAINERS_ALL"
#include "utils.hpp"
#include <numeric>

std::pair<std::size_t, std::size_t>
readFileSize(const std::string & fileName)
{
  std::pair<std::size_t, std::size_t> result;

  std::ifstream source;
  source.open(fileName, std::ios_base::in);
  std::string line;
  // read file line to get size of data in the file
  std::getline(source, line);
  std::istringstream in(line);
  std::string col1, col2;
  in >> col1; in >> col2;
  source.close();

  return std::make_pair(std::stoi(col1), std::stoi(col2));
}

template<typename sc_t, typename T>
std::vector<sc_t> computeScores(const T & rsv)
{
  const auto nR = rsv.rows();
  const auto k  = rsv.cols();

  std::vector<sc_t> result(nR);
  for (auto j=0; j<nR; ++j){
    for (auto xi=0; xi<k; ++xi){
      result[j] += rsv(j,xi)*rsv(j,xi);
    }
  }
  const auto scale = [&](sc_t & v) { v/=static_cast<sc_t>(k); };
  std::for_each(result.begin(), result.end(), scale);

  const auto print = [&](const sc_t & v) {std::cout << v << "\n"; };
  std::for_each(result.begin(), result.end(), print);

  // find min score
  const sc_t minValue = *std::min_element(result.begin(), result.end());

  sc_t sum = 0;
  auto add = [&](sc_t & v) { v+= std::abs(minValue); sum += v; };
  std::for_each(result.begin(), result.end(), add);

  auto normal = [&](sc_t & v) { v/=sum; };
  std::for_each(result.begin(), result.end(), normal);

  sum = std::accumulate(result.begin(), result.end(), 0.);
  std::cout << sum << std::endl;

  return result;
}

int main(int argc, char *argv[]){
  using sc_t		= double;
  using eig_dyn_vec	= Eigen::VectorXd;
  using dmat_t	        = Eigen::MatrixXd;

  if(argc==1){
    std::cout << "\nYou need to pass the folder where to read left singular vectors from.\n" << std::endl;
  }
  if(argc!=3){
    std::cout << "args: folder rank\n" << std::endl;
    return 1;
  }

  const auto dir = std::string(argv[1]);
  std::cout << "\nTarget Directory: " << dir << std::endl;
  const auto rank = std::atoi(argv[2]);
  std::cout << "\nTarget rank: " << rank << std::endl;

  // ---------------------------------------------------
  const std::vector<std::string> dofNames {"vp", "sp"};
  for (auto iDof=0; iDof<dofNames.size(); iDof++)
  {
    const auto thisDofName = dofNames[iDof];

    const std::string file = "test.txt";//dir + "/rsv_"+thisDofName+"_T.txt";
    std::cout << "\nReading snapshots: " << file << std::endl;

    const auto rsv = readDataIntoEigenMat<sc_t, int32_t, dmat_t>(file, rank, false);
    std::cout << "rsv size: " << rsv.rows() << " " << rsv.cols() << std::endl;

    const auto scores = computeScores<sc_t>(rsv);
  }
  return 0;
}
