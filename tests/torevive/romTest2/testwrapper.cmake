include(FindUnixCommands)

# remove possibly existing snapshots
execute_process(COMMAND ${BASH} -c "rm -rf snaps_vp snaps_sp lsv_vp lsv_sp")
execute_process(COMMAND ${BASH} -c "rm -rf final_generalized_coords_ascii*")
execute_process(COMMAND ${BASH} -c "rm -rf final_generalized_coords_bin*")

#--------------------------------------------
#--- first process things in ascii format ---
#--------------------------------------------
# 1. run the FOM exe
execute_process(COMMAND ${CMD_FOM} ${INPUT_FOM_ASCII} RESULT_VARIABLE CMD_RESULT)
message(${CMD_RESULT})
if(CMD_RESULT)
  message(FATAL_ERROR "Something wrong running FOM with ascii IO")
endif()

# 2. compute svd
execute_process(COMMAND ${CMD_SVD} ${CMAKE_CURRENT_SOURCE_DIR} 0 0 RESULT_VARIABLE CMD_RESULT)
message(${CMD_RESULT})
if(CMD_RESULT)
  message(FATAL_ERROR "Something wrong computing SVD with ascii IO")
endif()

# 3. run the ROM exe
execute_process(COMMAND ${CMD_ROM} ${INPUT_ROM_ASCII} RESULT_VARIABLE CMD_RESULT)
execute_process(COMMAND ${BASH} -c "cp final_generalized_coords_vp.txt final_generalized_coords_ascii.txt")
message(${CMD_RESULT})
if(CMD_RESULT)
  message(FATAL_ERROR "Something wrong running ROM with ascii IO")
endif()


#--------------------------------------------
#--- then process things in binary format ---
#--------------------------------------------
# 1. run the FOM exe
execute_process(COMMAND ${CMD_FOM} ${INPUT_FOM_BIN} RESULT_VARIABLE CMD_RESULT)
message(${CMD_RESULT})
if(CMD_RESULT)
  message(FATAL_ERROR "Something wrong running FOM with binary IO")
endif()

# 2. compute svd
execute_process(COMMAND ${CMD_SVD} ${CMAKE_CURRENT_SOURCE_DIR} 0 1 RESULT_VARIABLE CMD_RESULT)
message(${CMD_RESULT})
if(CMD_RESULT)
  message(FATAL_ERROR "Something wrong computing SVD with binary IO")
endif()

# 3. run the ROM exe
execute_process(COMMAND ${CMD_ROM} ${INPUT_ROM_BIN} RESULT_VARIABLE CMD_RESULT)
execute_process(COMMAND ${BASH} -c "cp final_generalized_coords_vp.txt final_generalized_coords_bin.txt")
message(${CMD_RESULT})
if(CMD_RESULT)
  message(FATAL_ERROR "Something wrong running ROM with binary IO")
endif()

#--------------------------------------------
#--- DO CHECKING
#--------------------------------------------
# check betwene ascii and binary
execute_process(COMMAND ${CMD_COMPARE}
  final_generalized_coords_bin.txt
  final_generalized_coords_ascii.txt
  1e-10
  RESULT_VARIABLE RES)
if(RES)
  message(FATAL_ERROR "Gen coords from binIO and asciiIO do not match")
endif()

# do final matching with gold values
execute_process(COMMAND ${CMD_COMPARE}
  final_generalized_coords_ascii.txt
  final_generalized_coords_gold.txt
  1e-10
  RESULT_VARIABLE RES)
if(RES)
  message(FATAL_ERROR "Gen coords from ascii do not match with gold")
endif()

# do final matching with gold values
execute_process(COMMAND ${CMD_COMPARE}
  final_generalized_coords_bin.txt
  final_generalized_coords_gold.txt
  1e-10
  RESULT_VARIABLE RES)
if(RES)
  message(FATAL_ERROR "Gen coords from bin do not match with gold")
endif()
