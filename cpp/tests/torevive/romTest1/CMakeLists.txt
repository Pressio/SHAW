
configure_file(input_fom1.yaml input_fom1.yaml COPYONLY)
configure_file(input_fom2.yaml input_fom2.yaml COPYONLY)

# configure_file(input.txt input.txt COPYONLY)
# configure_file(lsv_vp.txt lsv_vp.txt COPYONLY)
# configure_file(lsv_sp.txt lsv_sp.txt COPYONLY)
# configure_file(snaps_vp_gold.txt snaps_vp_gold.txt COPYONLY)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/../fullMesh101x201 DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# add_test(NAME romTest1EigenLTI
#   COMMAND ${CMAKE_COMMAND}
#   -DCMD_ROM=$<TARGET_FILE:shwave_eigen_rom>
#   -DCMD_COMPARE=$<TARGET_FILE:compareSnaps>
#   -DINPUT_FNAME=input.txt
#   -P ${CMAKE_CURRENT_SOURCE_DIR}/testwrapper.cmake
#   )

# if(HAVE_KOKKOS)
#   add_test(NAME romTest1KokkosLTI
#     COMMAND ${CMAKE_COMMAND}
#     -DCMD=$<TARGET_FILE:shwave_kokkos>
#     -DCMD2=$<TARGET_FILE:compareSnaps>
#     -DINPUT_FNAME=input.txt
#     -P ${CMAKE_CURRENT_SOURCE_DIR}/../testwrapper.cmake
#     )
# endif()
