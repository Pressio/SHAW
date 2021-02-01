
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../src/kokkos)

# test for input file with all sections present
configure_file(input1.yaml input1.yaml COPYONLY)
set(test_name parser_test_1)
add_executable(${test_name} main1.cc)
add_test(NAME ${test_name} COMMAND ${test_name} input1.yaml)
set_tests_properties(${test_name}
  PROPERTIES PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAILED"
  )

# test for input that does not have I/O section
configure_file(input2.yaml input2.yaml COPYONLY)
set(test_name parser_test_2)
add_executable(${test_name} main2.cc)
add_test(NAME ${test_name} COMMAND ${test_name} input2.yaml)
set_tests_properties(${test_name}
  PROPERTIES PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAILED"
  )
