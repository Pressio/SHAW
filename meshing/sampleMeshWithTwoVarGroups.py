#!/usr/bin/env python

import matplotlib.pyplot as plt
import sys, os, time
import numpy as np
from numpy import linspace, meshgrid
from matplotlib import cm
import collections

#from mesh_utils import *
#from unit_cell import UnitCell
#from earth_params import earthRadius, mantleThickness, cmbRadius

#--------------------------------------------------------------------------
def mainSampleMeshWithTwoVarGroups(workDir, targetPct,
                                   nth, nr, plotting,
                                   thL, thLDeg, thR, thRDeg,
                                   rSurf, rCmb, debugPrint,
                                   L, dth, dr):
  ## create list of cell objects for the full mesh
  # fullMeshCells = createCellsFullMesh(nth, nr, thL, thR, rSurf, rCmb, dth, dr)
  #if plotting != "none":
  #  for it in fullMeshCells: it.plotPoints(ax)

  # # store coordinates and label for each pt of the full grid
  # [coordsVp, coordsSp, dofLabelsSp] = storeCoordsAndLabelsFullMesh(fullMeshCells, nth, nr)

  # #   # total number of points in domain
  # totNumDofs = len(coordsVp) + len(coordsSp)

  # if (debugPrint == 1):
  #   print("----------------------------")
  #   print("--- graph for Vp ---")
  #   printDicPretty(gVpPrint)
  #   print("----------------------------")
  #   print("--- coeff for Vp ---")
  #   printDicPretty(coeffVpPrint)
  #   print("----------------------------")
  #   print("--- graph for Sp ---")
  #   printDicPretty(gSpPrint)

  # # -----------------------------------------------------
  # # number of pts where we compute residual or velocity
  # numResidualPtsVp = len(gVpPrint)
  # numResidualPtsSp = len(gSpPrint)
  # numResidualPts   = numResidualPtsVp + numResidualPtsSp
  # print ("numResidualPtsVp = ", numResidualPtsVp)
  # print ("numResidualPtsSp = ", numResidualPtsSp)
  # print ("numResidualPts = ", numResidualPts,
  #        " which is = ", numResidualPts/totNumDofs*100, " % of full mesh")

  #     if (targetPct < 0 ):
  #       print("For random sample mesh you need to set --pct=n")
  #       print("where n>0 and defines the % of full mesh to take")

  #     # get graphs for the sample mesh only
  #     [gVp, coeffVp, gSp] = buildGraphSampleMesh(fullMeshCells, nth, nr, targetPct)

  #     if (debugPrint == 1):
  #       print("----------------------------")
  #       print("--- graph for Vp ---")
  #       printDicPretty(gVp)
  #       print("----------------------------")
  #       print("--- coeff for Vp ---")
  #       printDicPretty(coeffVp)
  #       print("----------------------------")
  #       print("--- graph for Sp ---")
  #       printDicPretty(gSp)

  #     # Now we have graphs for the sample mesh cells for vp and srp/stp dofs.
  #     # These graphs contain the points where we need to compute the residual
  #     # as well as those where only the state is needed.
  #     # These graphs are defined in terms of the GIDs of the full mesh.
  #     # We need to uniquely enumerate all the sample mesh cells
  #     # separating the vp from srp/stp dofs.

  #     # We need to loop over the vp graph and the sp graph, and extract
  #     # all unique GIDs of the srp/stp nodes, since there might be duplicates
  #     # because of the connectivity.
  #     allVpGids = []
  #     allSpGids = []
  #     # loop over the vp graph first
  #     for k,v in gVp.items():
  #       # for vp, append key since it is a gid of a vp dof
  #       allVpGids.append(np.int64(k))
  #       # for Sp, append gids of the connected points since these in the vp graph
  #       # are either stp or srp dof.
  #       for j in v:
  #         allSpGids.append(np.int64(j))

  #     # loop over the graph for the srp/stp dofs
  #     for k,v in gSp.items():
  #       if (k>=0): allSpGids.append(np.int64(k))
  #       # dealing with srp dof, only take the north and south neighbors
  #       if (dofLabelsSp[k] == 1):
  #         allVpGids.append( np.int64(v[1]) )
  #         allVpGids.append( np.int64(v[3]) )
  #       # dealing with stp dof, only take the west and east neighbors
  #       if (dofLabelsSp[k] == 2):
  #         allVpGids.append( np.int64(v[0]) )
  #         allVpGids.append( np.int64(v[2]) )

  #     # remove duplicates and sort
  #     allVpGids = list(dict.fromkeys(allVpGids)); allVpGids.sort()
  #     allSpGids = list(dict.fromkeys(allSpGids)); allSpGids.sort()
  #     if (debugPrint == 1):
  #       print("----------------------------")
  #       print("Sample mesh unique GIDs for Vp:"); print(allVpGids)
  #       print("----------------------------")
  #       print("Sample mesh unique GIDs for srp/stp:"); print(allSpGids)
  #       print("\n")

  #     # -----------------------------------------------------
  #     # The sample mesh graph contains the GIDs
  #     # wrt the FULL mesh because we simply extracted a subset from the full mesh.
  #     # We enumerate the sample mesh points with new indexing
  #     # and create a map of full-mesh gids to new gids

  #     # fm_to_sm_map is such that:
  #     # - key   = GID wrt full mesh
  #     # - value = GID wrt sample mesh
  #     fm_to_sm_map_vp = collections.OrderedDict()
  #     fm_to_sm_map_vp = {pt:i for (pt,i) in zip(allVpGids, range(len(allVpGids)))}
  #     sm_to_fm_map_vp = {v:k for k,v in fm_to_sm_map_vp.items()}

  #     # sm_to_fm_map is such that:
  #     # - key   = GID wrt sample mesh
  #     # - value = GID wrt full mesh
  #     fm_to_sm_map_sp = collections.OrderedDict()
  #     fm_to_sm_map_sp = {pt:i for (pt,i) in zip(allSpGids, range(len(allSpGids)))}
  #     sm_to_fm_map_sp = {v:k for k,v in fm_to_sm_map_sp.items()}

  #     if (debugPrint == 1):
  #       print("----------------------------")
  #       print("--- fm_to_sm / sm_to_fm map for vp ---")
  #       printDicPretty2(fm_to_sm_map_vp, sm_to_fm_map_vp)
  #       print("----------------------------")
  #       print("--- fm_to_sm / sm_to_fm map for sp ---")
  #       printDicPretty2(fm_to_sm_map_sp, sm_to_fm_map_sp)

  #     # -----------------------------------------------------
  #     # print sm_to_fm file for Vp
  #     # -----------------------------------------------------
  #     f1 = open("sm_to_fm_vp.dat","w+")
  #     for k,v in sm_to_fm_map_vp.items():
  #       f1.write("%11d %11d\n" % (k, v))
  #     f1.close()

  #     # -----------------------------------------------------
  #     # print sm_to_fm file for Sp
  #     # -----------------------------------------------------
  #     f1 = open("sm_to_fm_sp.dat","w+")
  #     for k,v in sm_to_fm_map_sp.items():
  #       f1.write("%11d %11d\n" % (k, v))
  #     f1.close()

  #     # now that we have the fm <-> sm mapping, we can assemble the actual sample mesh
  #     # graph for vp that will be printed to file, which consists of the vp GIDs indexed wrt sample mesh
  #     print("------\n")
  #     gVpPrint     = collections.OrderedDict()
  #     coeffVpPrint = collections.OrderedDict()
  #     # loop over the original gVp graph (which still contains the GIDs wrt full mesh)
  #     for (k1,v1), (k2,v2) in zip(gVp.items(), coeffVp.items()):
  #       # k1 = k2 = vp dof global ID wrt full mesh enumeration
  #       # v1 = global IDs for the connecting dofs wrt full mesh enumeration
  #       # v2 = coeffcients for vp stencil at each point

  #       # use mapping to get the global ID of current vp dof wrt sample mesh indexing
  #       smGID = fm_to_sm_map_vp[k1]

  #       # loop over adjecent nodes (which are either srp or stp),
  #       # map their global ID to the ID wrt sample mesh
  #       # (remember the trick for ghost points which are stored as negative indices)
  #       # and store the converted index
  #       gVpPrint[smGID] = [ fm_to_sm_map_sp[i] for i in v1 ]

  #       # for the coeffs we keep everything as is since we only need to map the keys
  #       coeffVpPrint[smGID] = v2

  #     if (debugPrint == 1):
  #       print("----------------------------")
  #       printDicPretty2(gVpPrint, coeffVpPrint)


  #     # now that we have the fm_sm mapping and viceversa we can assemble the sample
  #     # graph for srp/stp that will be printed to file, which consists of GIDs wrt sample mesh
  #     print("------\n")
  #     gSpPrint = collections.OrderedDict()
  #     for k, v in gSp.items():
  #       # k = srp or stp dof global ID wrt full mesh enumeration
  #       # v = global IDs for the connecting dofs wrt full mesh enumeration

  #       # get the global ID wrt sample mesh indexing of current dof
  #       smGID = fm_to_sm_map_sp[k]

  #       stencilGIDs = v
  #       # if the dof being tested is srp, then we only need the gids of the
  #       # north and south vp nodes
  #       if (dofLabelsSp[k] == 1):
  #         stencilGIDs[1] = fm_to_sm_map_vp[v[1]]
  #         stencilGIDs[3] = fm_to_sm_map_vp[v[3]]

  #       # if the dof being tested is stp, then we only need the gids of the
  #       # west and east vp nodes
  #       if (dofLabelsSp[k] == 2):
  #         stencilGIDs[0] = fm_to_sm_map_vp[v[0]]
  #         stencilGIDs[2] = fm_to_sm_map_vp[v[2]]

  #       gSpPrint[smGID] = stencilGIDs


  #     if (debugPrint == 1):
  #       print("----------------------------")
  #       printDicPretty(gSpPrint)
  #       print("----------------------------")

  #     # -----------------------------------------------------
  #     # number of pts where we compute residual or velocity
  #     numResidualPtsVp = len(gVpPrint)
  #     numResidualPtsSp = len(gSpPrint)
  #     numResidualPts   = numResidualPtsVp + numResidualPtsSp
  #     print ("numResidualPtsVp = ", numResidualPtsVp)
  #     print ("numResidualPtsSp = ", numResidualPtsSp)
  #     print ("numResidualPts = ", numResidualPts,
  #            " which is = ", numResidualPts/totNumDofs*100, " % of full mesh")

  #     # -----------------------------------------------------
  #     # number of pts where we store the state
  #     numStatePtsVp = len(allVpGids)
  #     numStatePtsSp = len(allSpGids)
  #     numStatePts   = numStatePtsVp + numStatePtsSp
  #     print ("numStatePtsVp = ", numStatePtsVp)
  #     print ("numStatePtsSp = ", numStatePtsSp)
  #     print ("numStatePts = ", numStatePts,
  #            " which is = ", numStatePts/totNumDofs*100, " % of full mesh")

  #     # -----------------------------------------------------
  #     if plotting != "none" and samplingType == "random":
  #       th_fm = np.asarray([v[0] for (k,v) in coordsVp.items()])
  #       r_fm  = np.asarray([v[1] for (k,v) in coordsVp.items()])
  #       ax.scatter(th_fm, r_fm, marker='o', s=4, edgecolor='r', facecolor='r')
  #       th_fm = np.asarray([v[0] for (k,v) in coordsSp.items()])
  #       r_fm  = np.asarray([v[1] for (k,v) in coordsSp.items()])
  #       ax.scatter(th_fm, r_fm, marker='s', s=4, edgecolor='k', facecolor='k')

  #       # subselect coords for vp points in sample mesh
  #       th_fm = np.asarray([v[0] for (k,v) in coordsVp.items()])
  #       r_fm  = np.asarray([v[1] for (k,v) in coordsVp.items()])
  #       r_sm_vp, th_sm_vp  = r_fm[ list(fm_to_sm_map_vp.keys()) ], th_fm[ list(fm_to_sm_map_vp.keys()) ]
  #       ax1.scatter(th_sm_vp, r_sm_vp, marker='o', s=4, edgecolor='r', facecolor='r')
  #       # for i in range(len(allVpGids)):
  #       #   ax1.text(th_sm_vp[i], r_sm_vp[i], str( fm_to_sm_map_vp[allVpGids[i]] ),
  #       #          verticalalignment='top', horizontalalignment='center', fontsize=14, color='r')

  #       # subselect coords for srp/stp points in sample mesh
  #       th_fm = np.asarray([v[0] for (k,v) in coordsSp.items()])
  #       r_fm  = np.asarray([v[1] for (k,v) in coordsSp.items()])
  #       r_sm_sp, th_sm_sp  = r_fm[ list(fm_to_sm_map_sp.keys()) ], th_fm[ list(fm_to_sm_map_sp.keys()) ]
  #       ax1.scatter(th_sm_sp, r_sm_sp, marker='s', s=4, edgecolor='k', facecolor='k')
  #       # for i in range(len(allSpGids)):
  #       #   ax1.text(th_sm_sp[i], r_sm_sp[i], str( fm_to_sm_map_sp[allSpGids[i]] ),
  #       #          verticalalignment='top', horizontalalignment='center', fontsize=14, color='k')
