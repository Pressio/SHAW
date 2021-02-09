
configure_file(input.yaml input.yaml COPYONLY)

configure_file(snaps_vp_gold_0.txt snaps_vp_gold_0.txt COPYONLY)
configure_file(snaps_vp_gold_1.txt snaps_vp_gold_1.txt COPYONLY)
configure_file(snaps_vp_gold_2.txt snaps_vp_gold_2.txt COPYONLY)
configure_file(snaps_sp_gold_0.txt snaps_sp_gold_0.txt COPYONLY)
configure_file(snaps_sp_gold_1.txt snaps_sp_gold_1.txt COPYONLY)
configure_file(snaps_sp_gold_2.txt snaps_sp_gold_2.txt COPYONLY)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/../fullMesh21x51 DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_test(NAME multiPeriodsForcingRank1
  COMMAND ${CMAKE_COMMAND}
  -DCMD_FOM=$<TARGET_FILE:shwave_fom>
  -DCMD_COMPARE=$<TARGET_FILE:compareSnaps>
  -DINPUT_FNAME=input.yaml
  -P ${CMAKE_CURRENT_SOURCE_DIR}/test.cmake
  )
