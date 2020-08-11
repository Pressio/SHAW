
configure_file(input1.yaml input1.yaml COPYONLY)
configure_file(input2.yaml input2.yaml COPYONLY)

configure_file(snaps_vp_gold.txt snaps_vp_gold.txt COPYONLY)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/../fullMesh21x51 DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

if(HAVE_KOKKOS)
add_test(NAME fomSymmetryAxisThetaPiKokkos1
  COMMAND ${CMAKE_COMMAND}
  -DCMD_FOM=$<TARGET_FILE:shwave_fom>
  -DCMD_COMPARE=$<TARGET_FILE:compareSnaps>
  -DINPUT_FNAME=input1.yaml
  -P ${CMAKE_CURRENT_SOURCE_DIR}/../testwrapper_fom.cmake
  )

add_test(NAME fomSymmetryAxisThetaPiKokkos2
  COMMAND ${CMAKE_COMMAND}
  -DCMD_FOM=$<TARGET_FILE:shwave_fom>
  -DCMD_COMPARE=$<TARGET_FILE:compareSnaps>
  -DINPUT_FNAME=input2.yaml
  -P ${CMAKE_CURRENT_SOURCE_DIR}/../testwrapper_fom.cmake
  )
endif()
