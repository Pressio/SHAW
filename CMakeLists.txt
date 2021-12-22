
cmake_minimum_required(VERSION 3.16.0)
project(shaw CXX)

# guard against in-source builds
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds are not allowed.
Remove CMakeCache.txt. Make a separate build directory and run CMake from there.")
endif()

# default to release if build type is empty
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

# we need c++14
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#==================================

find_package(yaml-cpp REQUIRED)
message(STATUS ${YAML_CPP_INCLUDE_DIR})
include_directories(${YAML_CPP_INCLUDE_DIR})

find_package(KokkosKernels REQUIRED)

# executables
add_executable(
  shawExe
  ${CMAKE_CURRENT_SOURCE_DIR}/src/kokkos/main.cc)
target_link_libraries(shawExe dl ${YAML_CPP_LIBRARIES} Kokkos::kokkoskernels)

add_executable(
  extractStateFromSnaps
  ${CMAKE_CURRENT_SOURCE_DIR}/src/tools/main_extract_state_from_snaps.cc)
target_link_libraries(extractStateFromSnaps dl ${YAML_CPP_LIBRARIES} Kokkos::kokkoskernels)

add_executable(
  reconstructFomState
  ${CMAKE_CURRENT_SOURCE_DIR}/src/tools/main_reconstruct_fom_state.cc)
target_link_libraries(reconstructFomState dl ${YAML_CPP_LIBRARIES} Kokkos::kokkoskernels)

add_executable(
  reconstructSeismogram
  ${CMAKE_CURRENT_SOURCE_DIR}/src/tools/main_reconstruct_seismo.cc)
target_link_libraries(reconstructSeismogram dl ${YAML_CPP_LIBRARIES} Kokkos::kokkoskernels)

find_package(OpenMP)
add_executable(computeThinSVD ${CMAKE_CURRENT_SOURCE_DIR}/src/tools/main_eigen_svd.cc)
target_compile_options(computeThinSVD PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-march=native>)
target_include_directories(computeThinSVD PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/eigen)
target_link_libraries(computeThinSVD OpenMP::OpenMP_CXX)

# tests
enable_testing()
add_subdirectory(tests)


# REMOVE BELOW
#link_directories   (${BLAS_LIB_DIR} ${YAMLCPP_LIB_DIR})
#link_libraries	   (dl ${OMPLINK} yaml-cpp Kokkos::kokkoskernels ${BLASLIBNAME} gfortran)
#link_directories   (${YAMLCPP_LIB_DIR})
#link_libraries	   (dl yaml-cpp Kokkos::kokkoskernels) # Kokkos::BLAS Kokkos::LAPACK)
# set(exeName compareSnaps)
# add_executable(${exeName} ${CMAKE_CURRENT_SOURCE_DIR}/src/tools/main_cmp_snaps.cc)
# target_include_directories(${exeName} PRIVATE ${EIGEN_INCLUDE_DIR})
# target_compile_definitions(${exeName} PRIVATE SHW_ENABLE_TPL_EIGEN)
# target_compile_options	(${exeName} PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-march=native>)
# set(exeName computeThinSVD)
# add_executable		(${exeName} ${CMAKE_CURRENT_SOURCE_DIR}/src/tools/main_eigen_svd.cc)
# target_compile_options	(${exeName} PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-march=native>)
