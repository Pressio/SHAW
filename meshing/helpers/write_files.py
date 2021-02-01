
import numpy as np

########################
###### PRINT FILES #####
########################
# formatting float
ffmt = ".20f"

def printMeshInfoFileCommon(f, thLDeg, thRDeg, rCmb, rSurf, dth, dr):
  f.write(("thL   %"+ffmt+"\n") % thLDeg)
  f.write(("thR   %"+ffmt+"\n") % thRDeg)
  f.write(("rCmb  %"+ffmt+"\n") % rCmb)
  f.write(("rSurf %"+ffmt+"\n") % rSurf)
  f.write(("dth   %"+ffmt+"\n") % dth)
  f.write(("dr    %"+ffmt+"\n") % dr)

def printSampleMeshInfoFileTwoVarGroups(nth, nr,
                                        thLDeg, thRDeg, rCmb, rSurf, dth, dr,
                                        numStatePtsVp, numResidualPtsVp,
                                        numStatePtsSp, numResidualPtsSp):
  f = open("mesh_info.dat","w+")
  printMeshInfoFileCommon(f, thLDeg, thRDeg, rCmb, rSurf, dth, dr)
  # when we do sample mesh, we need to differentiate because
  # we have more points where we compute the state than those where
  # we compute the residual
  f.write(("numResidualPtsVp %d\n") % numResidualPtsVp)
  f.write(("numStatePtsVp    %d\n") % numStatePtsVp)
  f.write(("numResidualPtsSp %d\n") % numResidualPtsSp)
  f.write(("numStatePtsSp    %d\n") % numStatePtsSp)
  f.close()

def printFullMeshInfoFileTwoVarGroups(nth, nr,
                                      thLDeg, thRDeg, rCmb, rSurf, dth, dr,
                                      numStatePtsVp, numStatePtsSp):
  f = open("mesh_info.dat","w+")
  printMeshInfoFileCommon(f, thLDeg, thRDeg, rCmb, rSurf, dth, dr)
  f.write(("numPtsVp %d\n") % numStatePtsVp)
  f.write(("numPtsSp %d\n") % numStatePtsSp)
  f.write(("nth %d\n") % nth)
  f.write(("nr %d\n") % nr)
  f.close()


def printMeshInfoFileThreeVarGroups(samplingType, nth, nr,
                                    thLDeg, thRDeg, rCmb, rSurf, dth, dr,
                                    numStatePtsVp, numResidualPtsVp,
                                    numStatePtsSrp, numResidualPtsSrp,
                                    numStatePtsStp, numResidualPtsStp):
  f = open("mesh_info.dat","w+")
  printMeshInfoFileCommon(f, thLDeg, thRDeg, rCmb, rSurf, dth, dr)
  if (samplingType == "full"):
    # for full mesh, the numResidualPts* == numStatePts*
    # because there is no difference between residual poits and state points
    assert(numResidualPtsVp == numStatePtsVp)
    assert(numResidualPtsSrp == numStatePtsSrp)
    assert(numResidualPtsStp == numStatePtsStp)
    f.write(("numPtsVp %d\n") % numStatePtsVp)
    f.write(("numPtsSrp %d\n") % numStatePtsSrp)
    f.write(("numPtsStp %d\n") % numStatePtsStp)
    f.write(("nth %d\n") % nth)
    f.write(("nr %d\n") % nr)
  else:
    # when we do sample mesh, we need to differentiate because
    # we have more points where we compute the state than those where
    # we compute the residual
    f.write(("numResidualPtsVp %d\n") % numResidualPtsVp)
    f.write(("numStatePtsVp    %d\n") % numStatePtsVp)
    f.write(("numResidualPtsSrp %d\n") % numResidualPtsSrp)
    f.write(("numStatePtsSrp    %d\n") % numStatePtsSrp)
    f.write(("numResidualPtsStp %d\n") % numResidualPtsStp)
    f.write(("numStatePtsStp    %d\n") % numStatePtsStp)

  f.close()

def writeConnectivityFileFullMesh(G, coords, isOnSym, fileName):
  f = open(fileName,"w+")
  for k in sorted(G.keys()):
    # print the gid of this point
    f.write("%d " % k)
    # print if this point is on sym axis
    f.write("%d " % isOnSym[k])
    # print coordinates of this vp point
    f.write(("%"+ffmt+" ") % (coords[k][0]))
    f.write(("%"+ffmt+" ") % (coords[k][1]))
    # loop over the connectivity points (which are all of type either srp or stp)
    for i in G[k]: f.write("%d " % i)
    f.write("\n")
  f.close()

def writeVpStencilCoeffsFileFullMesh(coeffVpPrint, fileName):
  f = open(fileName,"w+")
  for k in sorted(coeffVpPrint.keys()):
    # print the gid of this point
    f.write("%d " % k)
    # loop over the connectivity points (which are all of type either srp or stp)
    for i in coeffVpPrint[k]:
      f.write("%.1f " % i)

    f.write("\n")
  f.close()

def writeConnectivityFileWithLabelsFullMesh(G, coords, labels, isOnSym, fileName):
  f = open(fileName,"w+")
  for k in sorted(G.keys()):
    # print the gid of this point
    f.write("%d " % k)
    # print the label of dof type of this point: 0=vp, 1=srp, 2=stp
    f.write("%1d " % labels[k])
    # print if this point is on sym axis
    f.write("%d " % isOnSym[k])

    # print coordinates of this vp point
    f.write(("%"+ffmt+" ") % (coords[k][0]))
    f.write(("%"+ffmt+" ") % (coords[k][1]))
    # loop over the connectivity points (which are all of type either srp or stp)
    for i in G[k]: f.write("%d " % i)
    f.write("\n")
  f.close()












#   # -----------------------------------------------------
#   # print connectivity file for Vp
#   # -----------------------------------------------------
#   f = open("graph_vp.dat","w+")
#   for k in sorted(gVpPrint.keys()):
#     # k = vp dof global ID wrt sample mesh (potentially)
#     # so get the global ID wrt full mesh indexing of current dof
#     if (samplingType=="full"):
#       kGID = k
#     elif (samplingType=="random"):
#       kGID = sm_to_fm_map_vp[k]

#     # print the gid of this point
#     f.write("%d " % k)
#     # print the label of dof type of this point: 0=vp
#     f.write("%1d " % 0)
#     # print coordinates of this vp point
#     f.write(("%"+ffmt+" ") % (coordsVp[kGID][0]))
#     f.write(("%"+ffmt+" ") % (coordsVp[kGID][1]))
#     # loop over the connectivity points (which are all of type either srp or stp)
#     for i in gVpPrint[k]:
#       f.write("%d " % i)
#       # pring theta and r for each connected point
#       if (samplingType=="full"):
#         iGID = i
#       elif (samplingType=="random"):
#         iGID = sm_to_fm_map_sp[i]

#       f.write(("%"+ffmt+" ") % (coordsSp[iGID][0]))
#       f.write(("%"+ffmt+" ") % (coordsSp[iGID][1]))

#     f.write("\n")
#   f.close()

#   # -----------------------------------------------------
#   # print coeffs file for Vp
#   # -----------------------------------------------------
#   f = open("coeff_vp.dat","w+")
#   for k in sorted(coeffVpPrint.keys()):
#     # k = vp dof global ID wrt sample mesh (potentially)
#     # so get the global ID wrt full mesh indexing of current dof
#     if (samplingType=="full"):
#       kGID = k
#     elif (samplingType=="random"):
#       kGID = sm_to_fm_map_vp[k]

#     # print the gid of this point
#     f.write("%d " % k)
#     # # print coordinates of this vp point
#     # f.write("%.15f " % (coordsVp[kGID][0]))
#     # f.write("%.15f " % (coordsVp[kGID][1]))
#     # loop over the connectivity points (which are all of type either srp or stp)
#     for i in coeffVpPrint[k]:
#       f.write("%.1f " % i)

#     f.write("\n")
#   f.close()

#   # -----------------------------------------------------
#   # print connectivity file for Sp
#   # -----------------------------------------------------
#   f = open("graph_sp.dat","w+")
#   for k in sorted(gSpPrint.keys()):
#     # k = sp dof global ID wrt sample mesh (potentially)
#     # so get the global ID wrt full mesh indexing of current dof
#     if (samplingType=="full"):
#       kGID = k
#     elif (samplingType=="random"):
#       kGID = sm_to_fm_map_sp[k]

#     # print the gid of this point
#     f.write("%d " % k)
#     # print the label of dof type of this point: 0=vp, 1=srp, 2=stp
#     f.write("%1d " % dofLabelsSp[kGID])
#     # print coordinates of this vp point
#     f.write(("%"+ffmt+" ") % (coordsSp[kGID][0]))
#     f.write(("%"+ffmt+" ") % (coordsSp[kGID][1]))
#     # loop over the connectivity points (which are all of type either vp)
#     for i in gSpPrint[k]:
#       f.write("%d " % i)

#       # pring theta and r for each connected point
#       if (samplingType=="full"):
#         iGID = i
#       elif (samplingType=="random"):
#         iGID = sm_to_fm_map_vp[i]
#       f.write(("%"+ffmt+" ") % (coordsVp[iGID][0]))
#       f.write(("%"+ffmt+" ") % (coordsVp[iGID][1]))
#     f.write("\n")
#   f.close()
#   # -----------------------------------------------------

#   # move all files to target directory
#   os.system("mv *.dat {}".format(workDir))

#   if plotting == "show":
#     plt.show()
#   elif plotting =="print":
#     figFM.savefig('fullFM.png')
#     if samplingType == "random":
#       figSM.savefig('fullSM.png')
