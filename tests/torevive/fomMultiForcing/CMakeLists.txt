
configure_file(input.yaml input.yaml COPYONLY)

configure_file(snaps_vp_gold.txt snaps_vp_gold.txt COPYONLY)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/../fullMesh50x200 DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

if(HAVE_KOKKOS)
add_test(NAME fomMultiForcing
  COMMAND ${CMAKE_COMMAND}
  -DCMD_FOM=$<TARGET_FILE:shawExe>
  -DCMD_COMPARE=$<TARGET_FILE:compareSnaps>
  -DINPUT_FNAME=input.yaml
  -P ${CMAKE_CURRENT_SOURCE_DIR}/../testwrapper_fom.cmake
  )
endif()
