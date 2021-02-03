
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../src/kokkos)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/../fullMesh21x51 DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

configure_file(input.yaml input.yaml COPYONLY)
set(test_name seismogram_test)
add_executable(${test_name} main.cc)
add_test(NAME ${test_name} COMMAND ${test_name} input.yaml)
set_tests_properties(${test_name}
  PROPERTIES PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAILED"
  )
