include(FindUnixCommands)

set(IDS 0 1 2 3 4 5 6 7 8)

# remove possibly existing snapshots
foreach(ID IN LISTS IDS)
 execute_process(COMMAND ${BASH} -c "rm -rf snaps_vp_${ID} snaps_sp_${ID} seismogram_${ID}")
endforeach()

# first run the exe
execute_process(COMMAND ${CMD_FOM} ${INPUT_FNAME} RESULT_VARIABLE CMD_RESULT)
message(${CMD_RESULT})
if(RES)
  message(FATAL_ERROR "Fom run failed")
endif()

set(FILES "snaps_vp;snaps_sp;seismogram")
foreach(FF IN LISTS FILES)
  foreach(RID IN LISTS IDS)
    set(tol 1e-13)
    if(${FF} MATCHES "snaps_sp")
      set(tol 1e-10)
    endif()

    set(finalArg 1)
    if(${FF} MATCHES "seismogram")
      set(finalArg 0)
    endif()

    set(CMD "python compare.py ${FF}_${RID} ${FF}_${RID}_gold ${tol} ${finalArg}")
    execute_process(COMMAND ${BASH} -c ${CMD} RESULT_VARIABLE RES)
    if(RES)
      message(FATAL_ERROR "Diff for ${FF} is not clean")
    endif()
  endforeach()
endforeach()

