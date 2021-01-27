#!/usr/bin/env python

import matplotlib.pyplot as plt
import sys, os, time
import numpy as np
from numpy import linspace, meshgrid
from matplotlib import cm
import collections
import random

sys.path.insert(0, './helpers')
from full_mesh_utils import *
from utils import *
from write_files import *
from build_graph_impl import *
from convert_graph_dic_to_sparse_matrix import *
from plot_utils import *

#--------------------------------------------------------------------------
def mainFullMeshWithTwoVarGroups(workDir, nth, nr, plotting,
                                 thL, thLDeg, thR, thRDeg,
                                 rSurf, rCmb, debugPrint,
                                 L, dth, dr):
  varGroups = 2

  # create list of cell objects for the full mesh
  fullMeshCells = createCellsFullMesh(nth, nr, thL, thR, rSurf, rCmb, dth, dr, varGroups)
  # if plotting != "none":
  #   for it in fullMeshCells: it.plotPoints(ax)

  # store coordinates and label for each pt of the full grid
  [coordsVp, coordsSp, dofLabelsSp, isOnSymVp, isOnSymSp] = storeCoordsAndLabelsFullMesh(fullMeshCells, nth, nr, varGroups)

  #   # total number of points in domain
  totNumDofs = len(coordsVp) + len(coordsSp)

  # get graphs for the sample mesh only
  targetCells = range(0, nth*nr)
  [gVp, coeffVp, gSp] = buildGraphImpl(fullMeshCells, nth, nr, targetCells, varGroups)

  if (debugPrint == 1):
    print("----------------------------")
    print("--- graph for Vp ---")
    printDicPretty(gVp)
    print("----------------------------")
    print("--- coeff for Vp ---")
    printDicPretty(coeffVp)
    print("----------------------------")
    print("--- graph for Sp ---")
    printDicPretty(gSp)

  # -----------------------------------------------------
  # number of pts where we store the state is the same as residual points
  numStatePtsVp = len(gVp)
  numStatePtsSp = len(gSp)
  numStatePts   = numStatePtsVp + numStatePtsSp
  print ("numPtsVp = ", numStatePtsVp)
  print ("numPtsSp = ", numStatePtsSp)
  print ("numPts = ", numStatePts)

  # -----------------------------------------------------
  # convert mesh graph to sparse matrix
  if (plotting != "none"):
    fig, (ax1, ax2) = plt.subplots(1,2);
    M = convertGraphDicToSparseMatrix(gVp)
    ax1.spy(M, aspect='auto',marker='s',markersize=5)
    ax1.set_aspect(1)
    ax1.set_title("gVp")

    M2 = convertGraphDicToSparseMatrix(gSp)
    ax2.spy(M2, aspect='auto',marker='s',markersize=5)
    ax2.set_aspect(1)
    ax2.set_title("gSp")

  printFullMeshInfoFileTwoVarGroups(nth, nr, thLDeg, thRDeg, rCmb, rSurf,
                                    dth, dr, numStatePtsVp, numStatePtsSp)

  writeConnectivityFileFullMesh(gVp, coordsVp, isOnSymVp, "graph_vp.dat")
  writeVpStencilCoeffsFileFullMesh(coeffVp, "coeff_vp.dat")
  writeConnectivityFileWithLabelsFullMesh(gSp, coordsSp,
                                          dofLabelsSp, isOnSymSp, "graph_sp.dat")


  # mesh files are copied to a subdir of the workdir
  myWorkDir = workDir + "/mesh" + str(nr) + "x" + str(nth)
  # check if workdir exits and if not create it
  if not os.path.exists(myWorkDir): os.system('mkdir -p ' + myWorkDir)

  # move all files to target directory
  os.system("mv *.dat {}".format(myWorkDir))

  # -----------------------------------------------------
  if plotting != "none":
    fig = plt.figure(0)
    ax = fig.add_subplot(111, projection='polar')
    ax.set_rlabel_position(180)
    ax.set_aspect(1.0)
    #plotEarthSurf(ax)
    #plotCMB(ax)

    # the minus sign for thetat below is needed to plot things so that
    # it shows the grid on the right half of the circle
    # TODO: fix labels of theta

    th_fm = -np.asarray([v[0] for (k,v) in coordsVp.items()])+np.pi*0.5
    r_fm  = np.asarray([v[1] for (k,v) in coordsVp.items()])
    gids  = np.asarray([k for (k,v) in coordsVp.items()])
    ax.scatter(th_fm, r_fm, marker='o', s=4, edgecolor='r', facecolor='r')
    for i in range(len(gids)):
      ax.text(th_fm[i], r_fm[i], str( gids[i] ),
              verticalalignment='top', horizontalalignment='center',
              fontsize=10, color='r')

    th_fm = -np.asarray([v[0] for (k,v) in coordsSp.items()])+np.pi*0.5
    r_fm  = np.asarray([v[1] for (k,v) in coordsSp.items()])
    gids  = np.asarray([k for (k,v) in coordsSp.items()])
    ax.scatter(th_fm, r_fm, marker='s', s=4, edgecolor='k', facecolor='k')
    for i in range(len(gids)):
      ax.text(th_fm[i], r_fm[i], str( gids[i] ),
              verticalalignment='top', horizontalalignment='center',
              fontsize=10, color='k')
  plt.show()
