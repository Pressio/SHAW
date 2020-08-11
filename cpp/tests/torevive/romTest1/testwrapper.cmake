include(FindUnixCommands)

# remove possibly existing snapshots
execute_process(COMMAND ${BASH} -c "rm -rf snaps_vp snaps_sp")

# first run the exe
execute_process(COMMAND ${CMD_ROM} ${INPUT_FNAME} RESULT_VARIABLE CMD_RESULT)
message(${CMD_RESULT})
if(RES)
  message(FATAL_ERROR "run failed")
endif()

execute_process(COMMAND ${CMD_COMPARE} snaps_vp snaps_vp_gold.txt 1e-12 RESULT_VARIABLE RES)
message(${RES})
if(RES)
  message(FATAL_ERROR "Diff for vp is not clean")
endif()

