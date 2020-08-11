include(FindUnixCommands)

# remove possibly existing snapshots
execute_process(COMMAND ${BASH} -c "rm -rf snaps_vp snaps_sp")

# first run the exe
execute_process(COMMAND ${CMD_FOM} ${INPUT_FNAME} RESULT_VARIABLE CMD_RESULT)
message(${CMD_RESULT})
if(RES)
  message(FATAL_ERROR "Fom run failed")
endif()

execute_process(COMMAND ${CMD_COMPARE} snaps_vp_0 snaps_vp_gold.txt 1e-13 RESULT_VARIABLE RES)
message(${RES})
if(RES)
  message(FATAL_ERROR "Diff for vp is not clean")
endif()




#execute_process(COMMAND ${BASH} -c "diff -b snaps_vp.txt snaps_vp_gold.txt" RESULT_VARIABLE RES)
# execute_process(COMMAND cmp --silent snaps_vp.txt snaps_vp_gold.txt || echo "PASSED")
# macro(EXEC_CHECK CMD)
# # if(CMD_RESULT)
# #     message(FATAL_ERROR "Error running ${CMD}")
# # else()
# #     if (BASH)
# #         execute_process(COMMAND ${BASH} -c "diff -b ${TEST_DATA_DIR}/test_input.txt ${TEST_DATA_DIR}/test_output.txt" RESULT_VARIABLE RES)
# #         if(RES)
# #             message(FATAL_ERROR "Diff is not clean")
# #         endif()
# #     else(BASH)
# #         message(FATAL_ERROR "BASH not found : no diff script run")
# #     endif(BASH)
# # endif()
# endmacro()
# exec_check(${CMD1})
