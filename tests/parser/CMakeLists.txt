
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

# test for bilayer model
configure_file(input3.yaml input3.yaml COPYONLY)
set(test_name parser_test_3)
add_executable(${test_name} main3.cc)
add_test(NAME ${test_name} COMMAND ${test_name} input3.yaml)
set_tests_properties(${test_name}
  PROPERTIES PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAILED"
  )

# test for prem
configure_file(input4.yaml input4.yaml COPYONLY)
set(test_name parser_test_4)
add_executable(${test_name} main4.cc)
add_test(NAME ${test_name} COMMAND ${test_name} input4.yaml)
set_tests_properties(${test_name}
  PROPERTIES PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAILED"
  )

# test for multiforcing, rank-1, sampling depth
configure_file(input5.yaml input5.yaml COPYONLY)
set(test_name parser_test_5)
add_executable(${test_name} main5.cc)
add_test(NAME ${test_name} COMMAND ${test_name} input5.yaml)
set_tests_properties(${test_name}
  PROPERTIES PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAILED"
  )

# test for multiforcing, rank-1, sampling periods
configure_file(input6.yaml input6.yaml COPYONLY)
set(test_name parser_test_6)
add_executable(${test_name} main6.cc)
add_test(NAME ${test_name} COMMAND ${test_name} input6.yaml)
set_tests_properties(${test_name}
  PROPERTIES PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAILED"
  )
