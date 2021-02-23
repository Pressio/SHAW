
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../compare.py compare.py COPYONLY)

configure_file(input.yaml input.yaml COPYONLY)

set(IDS 0 1 2)
foreach(ID IN LISTS IDS)
	configure_file(seismogram_${ID}_gold seismogram_${ID}_gold COPYONLY)
	configure_file(snaps_vp_${ID}_gold snaps_vp_${ID}_gold COPYONLY)
	configure_file(snaps_sp_${ID}_gold snaps_sp_${ID}_gold COPYONLY)
endforeach()

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/../fullMesh21x51 DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_test(NAME multiDepthsForcingRank1
  COMMAND ${CMAKE_COMMAND}
  -DCMD_FOM=$<TARGET_FILE:shwave_fom>
  -DINPUT_FNAME=input.yaml
  -P ${CMAKE_CURRENT_SOURCE_DIR}/test.cmake
  )
