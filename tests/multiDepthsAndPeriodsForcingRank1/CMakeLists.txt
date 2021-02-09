
configure_file(input.yaml input.yaml COPYONLY)

configure_file(snaps_vp_gold_0 snaps_vp_gold_0 COPYONLY)
configure_file(snaps_vp_gold_1 snaps_vp_gold_1 COPYONLY)
configure_file(snaps_vp_gold_2 snaps_vp_gold_2 COPYONLY)
configure_file(snaps_vp_gold_3 snaps_vp_gold_3 COPYONLY)
configure_file(snaps_vp_gold_4 snaps_vp_gold_4 COPYONLY)
configure_file(snaps_vp_gold_5 snaps_vp_gold_5 COPYONLY)
configure_file(snaps_vp_gold_6 snaps_vp_gold_6 COPYONLY)
configure_file(snaps_vp_gold_7 snaps_vp_gold_7 COPYONLY)
configure_file(snaps_vp_gold_8 snaps_vp_gold_8 COPYONLY)

configure_file(snaps_sp_gold_0 snaps_sp_gold_0 COPYONLY)
configure_file(snaps_sp_gold_1 snaps_sp_gold_1 COPYONLY)
configure_file(snaps_sp_gold_2 snaps_sp_gold_2 COPYONLY)
configure_file(snaps_sp_gold_3 snaps_sp_gold_3 COPYONLY)
configure_file(snaps_sp_gold_4 snaps_sp_gold_4 COPYONLY)
configure_file(snaps_sp_gold_5 snaps_sp_gold_5 COPYONLY)
configure_file(snaps_sp_gold_6 snaps_sp_gold_6 COPYONLY)
configure_file(snaps_sp_gold_7 snaps_sp_gold_7 COPYONLY)
configure_file(snaps_sp_gold_8 snaps_sp_gold_8 COPYONLY)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/../fullMesh21x51 DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_test(NAME multiPeriodsForcingRank1
  COMMAND ${CMAKE_COMMAND}
  -DCMD_FOM=$<TARGET_FILE:shwave_fom>
  -DCMD_COMPARE=$<TARGET_FILE:compareSnaps>
  -DINPUT_FNAME=input.yaml
  -P ${CMAKE_CURRENT_SOURCE_DIR}/test.cmake
  )
