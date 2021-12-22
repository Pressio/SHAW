
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

find_package(yaml-cpp REQUIRED HINTS ${YAMLCPP_DIR})
include_directories(${YAMLCPP_DIR}/include)

if(EXISTS "${YAMLCPP_DIR}/lib64/libyaml-cpp.a")
  set(YAML_CPP_LIBRARIES ${YAMLCPP_DIR}/lib64/libyaml-cpp.a)
endif()
if(EXISTS "${YAMLCPP_DIR}/lib/libyaml-cpp.a")
  set(YAML_CPP_LIBRARIES ${YAMLCPP_DIR}/lib/libyaml-cpp.a)

find_package(KokkosKernels REQUIRED HINTS ${KOKKOSKERNELS_DIR})

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
