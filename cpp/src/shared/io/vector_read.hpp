
#ifndef VECTOR_READ_HPP_
#define VECTOR_READ_HPP_

#include <fstream>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <vector>

template <typename T>
typename std::enable_if< is_kokkos_1dview<T>::value >::type
readAsciiVectorWithSize(const std::string fileName, T v)
{
  static_assert( is_accessible_on_host<T>::value,
		 "readAsciiVectorWithSize: the kokkos view must have HostSpace to read");

  std::ifstream source; source.open(fileName, std::ios_base::in);
  std::string line, colv;
  {
    std::getline(source, line);
    std::istringstream in(line);
    std::string col1;
    in >> col1;
    Kokkos::resize(v, std::stoi(col1));
  }

  // then read the actual data
  std::size_t iRow = 0;
  while (std::getline(source, line) ){
    std::istringstream in(line);
    in >> colv;
    v(iRow) = atof(colv.c_str());
    iRow++;
  }
  source.close();
}

#endif
