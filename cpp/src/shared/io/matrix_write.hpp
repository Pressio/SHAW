
#ifndef MATRIX_WRITE_HPP_
#define MATRIX_WRITE_HPP_

#include <fstream>
#include <iomanip>
#include <iostream>
#include <ostream>

namespace impl
{

template<typename sc_t>
void write_contig_matrix_to_binary(const std::string filename,
			    const sc_t * A,
			    std::size_t m,
			    std::size_t n,
			    bool printSize = true)
{
  std::ofstream out(filename, std::ios::out | std::ios::binary | std::ios::trunc);
  if (printSize){
    out.write((char*) (&m), sizeof(std::size_t));
    out.write((char*) (&n), sizeof(std::size_t));
  }
  out.write((char*) A, m*n*sizeof(sc_t) );
  out.close();
}

template<typename sc_t>
void write_contig_3d_view_to_binary(const std::string filename,
				    const sc_t * A,
				    std::size_t s0,
				    std::size_t s1,
				    std::size_t s2,
				    bool printSize = true)
{
  std::ofstream out(filename, std::ios::out | std::ios::binary | std::ios::trunc);
  if (printSize){
    out.write((char*) (&s0), sizeof(std::size_t));
    out.write((char*) (&s1), sizeof(std::size_t));
    out.write((char*) (&s2), sizeof(std::size_t));
  }
  out.write((char*) A, s0*s1*s2*sizeof(sc_t) );
  out.close();
}


template<typename mat_t>
void write_contig_3d_view_to_ascii(const std::string fileName,
				    const mat_t & A,
				    std::size_t s0,
				    std::size_t s1,
				    std::size_t s2,
				    bool printSize = true)
{
  std::ofstream file; file.open(fileName);
  if (printSize)
    file << s0 << " " << s1 << " " << s2 << std::endl;

  for (std::size_t k=0; k<s2; k++){
    for (std::size_t i=0; i<s0; i++){
      for (std::size_t j=0; j<s1; j++){
	file << std::setprecision(dblFmt) << A(i,j,k) << " ";
      }
      file << std::endl;
    }
  }
  file.close();
}

template <typename mat_t>
void write_matrix_to_ascii(const std::string fileName,
			   const mat_t & A,
			   std::size_t m,
			   std::size_t n,
			   bool printSize = true)
{
  std::ofstream file; file.open(fileName);
  if (printSize)
    file << m << " " << n << std::endl;

  for (std::size_t i=0; i<m; i++){
    for (std::size_t j=0; j<n; j++){
      file << std::setprecision(dblFmt) << A(i,j) << " ";
    }
    file << std::endl;
  }
  file.close();
}

}// impl namespace


// *** EIGEN  ***
template<class T>
typename std::enable_if<is_dynamic_matrix_eigen<T>::value>::type
writeToFile(const std::string fileName,
	    const T & A,
	    const bool useBinary,
	    const bool printSize = true)
{
  if (useBinary == 1)
    impl::write_contig_matrix_to_binary(fileName, A.data(), A.rows(), A.cols(), printSize);
  else
    impl::write_matrix_to_ascii(fileName, A, A.rows(), A.cols(), printSize);
}

// *** KOKKOS  ***
template<class T>
typename std::enable_if<is_kokkos_2dview<T>::value>::type
writeToFile(const std::string fileName,
	    const T A,
	    const bool useBinary,
	    const bool printSize = true)
{
  static_assert(is_accessible_on_host<T>::value,
		"cannot call writeToFile for view not accessible on host");

  if (useBinary == 1){
    if (!A.span_is_contiguous())
      throw std::runtime_error("Cannot write to binary a non-contig 2d view");

    impl::write_contig_matrix_to_binary(fileName, A.data(), A.extent(0), A.extent(1), printSize);
  }
  else
    impl::write_matrix_to_ascii(fileName, A, A.extent(0), A.extent(1), printSize);
}

template<class T>
typename std::enable_if<is_kokkos_3dview<T>::value>::type
writeToFile(const std::string fileName,
	    const T A,
	    const bool useBinary,
	    const bool printSize = true)
{
  static_assert(is_accessible_on_host<T>::value,
		"cannot call writeToFile for view not accessible on host");

  if (!A.span_is_contiguous())
    throw std::runtime_error("Cannot write to binary a non-contig 3d view");

  if (useBinary == 1){
    impl::write_contig_3d_view_to_binary(fileName, A.data(),
					 A.extent(0), A.extent(1), A.extent(2), printSize);
  }
  else
    impl::write_contig_3d_view_to_ascii(fileName, A,
					 A.extent(0), A.extent(1), A.extent(2), printSize);
}

#endif
