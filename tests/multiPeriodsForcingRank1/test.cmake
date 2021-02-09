include(FindUnixCommands)

# remove possibly existing snapshots
execute_process(COMMAND ${BASH} -c "rm -rf snaps_vp_0 snaps_sp_0")
execute_process(COMMAND ${BASH} -c "rm -rf snaps_vp_1 snaps_sp_1")
execute_process(COMMAND ${BASH} -c "rm -rf snaps_vp_2 snaps_sp_2")

# first run the exe
execute_process(COMMAND ${CMD_FOM} ${INPUT_FNAME} RESULT_VARIABLE CMD_RESULT)
message(${CMD_RESULT})
if(RES)
  message(FATAL_ERROR "Fom run failed")
endif()

execute_process(COMMAND ${CMD_COMPARE} snaps_vp_0 snaps_vp_gold_0.txt 1e-13 RESULT_VARIABLE RES)
if(RES)
  message(FATAL_ERROR "Diff for vp_0 is not clean")
endif()
execute_process(COMMAND ${CMD_COMPARE} snaps_sp_0 snaps_sp_gold_0.txt 1e-13 RESULT_VARIABLE RES)
if(RES)
  message(FATAL_ERROR "Diff for sp_0 is not clean")
endif()


execute_process(COMMAND ${CMD_COMPARE} snaps_vp_1 snaps_vp_gold_1.txt 1e-13 RESULT_VARIABLE RES)
if(RES)
  message(FATAL_ERROR "Diff for vp_1 is not clean")
endif()
execute_process(COMMAND ${CMD_COMPARE} snaps_sp_1 snaps_sp_gold_1.txt 1e-13 RESULT_VARIABLE RES)
if(RES)
  message(FATAL_ERROR "Diff for sp_1 is not clean")
endif()


execute_process(COMMAND ${CMD_COMPARE} snaps_vp_2 snaps_vp_gold_2.txt 1e-13 RESULT_VARIABLE RES)
if(RES)
  message(FATAL_ERROR "Diff for vp_2 is not clean")
endif()
execute_process(COMMAND ${CMD_COMPARE} snaps_sp_2 snaps_sp_gold_2.txt 1e-13 RESULT_VARIABLE RES)
if(RES)
  message(FATAL_ERROR "Diff for sp_2 is not clean")
endif()
