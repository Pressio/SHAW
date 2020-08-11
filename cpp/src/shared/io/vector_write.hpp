
#ifndef VECTOR_WRITE_HPP_
#define VECTOR_WRITE_HPP_

#include <fstream>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <vector>

namespace impl
{

template<typename sc_t, typename size_t>
void write_vector_to_binary(const std::string filename,
			    const sc_t * A,
			    size_t n,
			    bool printSize = true)
{
  std::ofstream out(filename, std::ios::out | std::ios::binary | std::ios::trunc);
  if (printSize)
    out.write((char*) (&n), sizeof(size_t));

  out.write((char*) A, n*sizeof(sc_t) );
  out.close();
}

template <typename mat_t, typename size_t>
void write_vector_to_ascii(const std::string fileName,
			   const mat_t & A,
			   size_t n,
			   bool printSize = true)
{
  std::ofstream file; file.open(fileName);
  if (printSize)
    file << n << std::endl;

  for (size_t i=0; i<n; i++){
    file << std::setprecision(dblFmt) << A(i) << " \n";
  }
  file.close();
}

}// impl namespace

template <typename T>
typename std::enable_if< is_kokkos_1dview<T>::value >::type
writeToFile(const std::string fileName,
	    const T & v,
	    const bool useBinary,
	    bool printSize = true)
{
  static_assert(is_accessible_on_host<T>::value,
		"cannot call writeToFile for view not accessible on host");

  if (useBinary == 1)
    impl::write_vector_to_binary(fileName, v.data(), v.extent(0), printSize);
  else
    impl::write_vector_to_ascii(fileName, v, v.extent(0), printSize);
}

#endif
