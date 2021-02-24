
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../compare.py compare.py COPYONLY)

configure_file(input.yaml input.yaml COPYONLY)

configure_file(snaps_vp_0_gold snaps_vp_0_gold COPYONLY)
configure_file(snaps_sp_0_gold snaps_sp_0_gold COPYONLY)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/../fullMesh21x51 DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_test(NAME fomNearCmb
  COMMAND ${CMAKE_COMMAND}
  -DCMD_FOM=$<TARGET_FILE:shawExe>
  -DINPUT_FNAME=input.yaml
  -P ${CMAKE_CURRENT_SOURCE_DIR}/test.cmake
  )
