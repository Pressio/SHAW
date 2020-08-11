
configure_file(input_fom_ascii.txt input_fom_ascii.txt COPYONLY)
configure_file(input_rom_ascii.txt input_rom_ascii.txt COPYONLY)
configure_file(input_fom_bin.txt input_fom_bin.txt COPYONLY)
configure_file(input_rom_bin.txt input_rom_bin.txt COPYONLY)
configure_file(final_generalized_coords_gold.txt final_generalized_coords_gold.txt COPYONLY)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/../fullMesh101x201 DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_test(NAME romTest1EigenLTI
  COMMAND ${CMAKE_COMMAND}
  -DCMD_FOM=$<TARGET_FILE:shwave_eigen_fom>
  -DCMD_ROM=$<TARGET_FILE:shwave_eigen_rom>
  -DCMD_SVD=$<TARGET_FILE:computeThinSVD>
  -DCMD_COMPARE=$<TARGET_FILE:compareSnaps>
  -DINPUT_FOM_ASCII=input_fom_ascii.txt
  -DINPUT_ROM_ASCII=input_rom_ascii.txt
  -DINPUT_FOM_BIN=input_fom_bin.txt
  -DINPUT_ROM_BIN=input_rom_bin.txt
  -P ${CMAKE_CURRENT_SOURCE_DIR}/testwrapper.cmake
  )

# if(HAVE_KOKKOS)
#   add_test(NAME romTest1KokkosLTI
#   COMMAND ${CMAKE_COMMAND}
#   -DCMD=$<TARGET_FILE:shwave_kokkos>
#   -DCMD2=$<TARGET_FILE:computeThinSVD>
#   -DCMD3=$<TARGET_FILE:compareSnaps>
#   -DINPUT_FOM_ASCII=input_fom_ascii.txt
#   -DINPUT_ROM_ASCII=input_rom_ascii.txt
#   -DINPUT_FOM_BIN=input_fom_bin.txt
#   -DINPUT_ROM_BIN=input_rom_bin.txt
#   -P ${CMAKE_CURRENT_SOURCE_DIR}/testwrapper.cmake
#   )
# endif()
