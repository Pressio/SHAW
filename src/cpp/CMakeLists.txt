
cmake_minimum_required(VERSION 3.11.0)
project(shwave CXX)

# guard against in-source builds
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds are not allowed.
Make a separate build directory and run CMake from there. Remove CMakeCache.txt.")
endif()

# default to release if build type is empty
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

# we need c++14
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# if(${HAVE_OMP})
#   message("You enabled OpenMP")
#   find_package(OpenMP)
#   set(OMPLINK OpenMP::OpenMP_CXX)
# else()
#   set(OMPLINK "")
# endif()

#==================================
# for every target we need:
find_package	   (KokkosKernels REQUIRED)
include_directories(${EIGEN_INCLUDE_DIR} ${YAMLCPP_INCLUDE_DIR} ${KOKKOS_INCLUDE_DIR})
link_directories   (${BLAS_LIB_DIR} ${YAMLCPP_LIB_DIR})
#link_libraries	   (dl ${OMPLINK} yaml-cpp Kokkos::kokkoskernels ${BLASLIBNAME} gfortran)
link_libraries	   (dl yaml-cpp Kokkos::kokkoskernels) # Kokkos::BLAS Kokkos::LAPACK)

#==================================

set(exeName	shwave_fom)
add_executable	(${exeName} ${CMAKE_CURRENT_SOURCE_DIR}/src/kokkos/main_fom.cc)

# set(exeName	shwave_rom)
# add_executable  (${exeName} ${CMAKE_CURRENT_SOURCE_DIR}/src/kokkos/main_rom.cc)

# set(exeName reconstructFomState)
# add_executable	(${exeName} ${CMAKE_CURRENT_SOURCE_DIR}/src/tools/main_reconstruct_fom_state.cc)

# set(exeName reconstructSeismogram)
# add_executable	(${exeName} ${CMAKE_CURRENT_SOURCE_DIR}/src/tools/main_reconstruct_seismo.cc)

# set(exeName extractStateFromSnaps)
# add_executable	(${exeName} ${CMAKE_CURRENT_SOURCE_DIR}/src/tools/main_extract_state_from_snaps.cc)

set(exeName compareSnaps)
add_executable		(${exeName} ${CMAKE_CURRENT_SOURCE_DIR}/src/tools/main_cmp_snaps.cc)
target_compile_options	(${exeName} PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-march=native>)

set(exeName computeThinSVD)
add_executable		(${exeName} ${CMAKE_CURRENT_SOURCE_DIR}/src/tools/main_eigen_svd.cc)
target_compile_options	(${exeName} PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-march=native>)

#==================================
# tell CMake to use CTest extension
enable_testing()
# add subdir where the unit tests are
add_subdirectory(tests)
#==================================


# set(exeName gemm)
# add_executable (${exeName} ${CMAKE_CURRENT_SOURCE_DIR}/src/benchmark/gemm.cc)
# #target_compile_definitions(${exeName} PUBLIC -DEIGEN_USE_BLAS)
# set(exeName gemv)
# add_executable (${exeName} ${CMAKE_CURRENT_SOURCE_DIR}/src/benchmark/gemv.cc)
# set(exeName spmv)
# add_executable (${exeName} ${CMAKE_CURRENT_SOURCE_DIR}/src/benchmark/main_spmv.cc)
