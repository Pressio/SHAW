
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/../fullMesh4x5 DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

configure_file(input.yaml input.yaml COPYONLY)
set(test_name stress_labels)
add_executable(${test_name} main.cc)

add_test(NAME ${test_name} COMMAND ${test_name} input.yaml)
set_tests_properties(${test_name}
  PROPERTIES PASS_REGULAR_EXPRESSION "PASS"
  FAIL_REGULAR_EXPRESSION "FAILED"
  )
