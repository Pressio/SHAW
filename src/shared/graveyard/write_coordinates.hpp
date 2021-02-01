
#ifndef WRITE_COORDINATES_HPP_
#define WRITE_COORDINATES_HPP_

// #include "vector_io.hpp"

namespace
{
// template <typename sc_t, typename vec_t>
// void writeCoordinatesToFile(const vec_t & thVp,
// 			    const vec_t & rrVpInv,
// 			    const std::size_t extent1,
// 			    const vec_t & thSp,
// 			    const vec_t & rrSpInv,
// 			    const std::size_t extent2)
// {
//   constexpr auto one = constants<sc_t>::one();
//   std::ofstream file;
//   file.open("coords_vp.txt");
//   for(auto i=0; i < extent1; i++){
//     file << std::setprecision(dblFmt)
// 	 << thVp[i] << " " << one/rrVpInv[i] << std::endl;
//   }
//   file.close();

//   file.open("coords_sp.txt");
//   for(auto i=0; i < extent2; i++){
//     file << std::setprecision(dblFmt)
// 	 << thSp[i] << " " << one/rrSpInv[i] << std::endl;
//   }
//   file.close();
// }

template <typename sc_t, typename coords_t>
void writeCoordinatesToFile(const std::string fileName,
			    const coords_t & coords,
			    const std::size_t extent1)
{
  constexpr auto one = constants<sc_t>::one();
  std::ofstream file; file.open(fileName);
  for(auto i=0; i < extent1; i++){
    file << std::setprecision(dblFmt)
	 << coords(i,0) << " " << one/coords(i,1)
	 << std::endl;
  }
  file.close();
}

}//anonym namespace

// template <typename vec_t>
// typename std::enable_if<is_vector_eigen<vec_t>::value>::type
// writeCoordinatesToFile(const vec_t & thVp,
// 		       const vec_t & rrVpInv,
// 		       const vec_t & thSp,
// 		       const vec_t & rrSpInv)
// {
//   using sc_t = typename vec_t::Scalar;
//   writeCoordinatesToFile<sc_t>(thVp, rrVpInv, thVp.size(),
// 			       thSp, rrSpInv, thSp.size());
// }

template <typename coords_t>
typename std::enable_if<is_kokkos_view<coords_t>::value>::type
writeCoordinatesToFile(const coords_t & coords, const dofId dof)
{
  using sc_t = typename coords_t::traits::value_type;
  const std::string dofName = dofIdToString(dof);
  const std::string filePath = "coords_" + dofName + ".txt";
  writeCoordinatesToFile<sc_t>(filePath, coords, coords.extent(0));
}

#endif
