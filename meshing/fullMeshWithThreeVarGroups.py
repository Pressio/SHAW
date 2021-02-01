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
def mainFullMeshWithThreeVarGroups(workDir, nth, nr, plotting,
                                   thL, thLDeg, thR, thRDeg,
                                   rSurf, rCmb, debugPrint,
                                   L, dth, dr):
  varGroups = 3

  # create list of cell objects for the full mesh
  fullMeshCells = createCellsFullMesh(nth, nr, thL, thR, rSurf, rCmb, dth, dr, varGroups)
  # #if plotting != "none":
  # #  for it in fullMeshCells: it.plotPoints(ax)

  # store coordinates and label for each pt of the full grid
  [coordsVp, coordsSrp, coordsStp] = storeCoordsAndLabelsFullMesh(fullMeshCells, nth, nr, varGroups)

  # total number of points in domain
  totNumDofs = len(coordsVp) + len(coordsSrp) + len(coordsStp)

  # get graphs for the sample mesh only
  targetCells = range(0, nth*nr)
  [gVpSrp, gVpStp, coeffVpSrp, coeffVpStp, gSrp, gStp] = buildGraphImpl(fullMeshCells, nth, nr, targetCells, varGroups)
  assert(len(gVpSrp) == len(gVpStp))

  if (debugPrint == 1):
    print("----------------------------")
    print("--- graph for Vp wrt srp ---")
    printDicPretty(gVpSrp)
    print("----------------------------")
    print("--- graph for Vp wrt stp ---")
    printDicPretty(gVpStp)
    print("----------------------------")
    print("--- coeff for Vp wrt srp ---")
    printDicPretty(coeffVpSrp)
    print("----------------------------")
    print("--- coeff for Vp wrt stp ---")
    printDicPretty(coeffVpStp)
    print("----------------------------")
    print("--- graph for Srp ---")
    printDicPretty(gSrp)
    print("----------------------------")
    print("--- graph for Stp ---")
    printDicPretty(gStp)

  # -----------------------------------------------------
  # number of pts where we compute residual or velocity
  numResidualPtsVp = len(gVpSrp)
  numResidualPtsSrp = len(gSrp)
  numResidualPtsStp = len(gStp)
  numResidualPts   = numResidualPtsVp + numResidualPtsSrp + numResidualPtsStp
  print ("numResidualPtsVp = ", numResidualPtsVp)
  print ("numResidualPtsSrp = ", numResidualPtsSrp)
  print ("numResidualPtsStp = ", numResidualPtsStp)
  print ("numResidualPts = ", numResidualPts,
         " which is = ", numResidualPts/totNumDofs*100, " % of full mesh")

  # -----------------------------------------------------
  # convert mesh graph to sparse matrix
  if (plotting != "none"):
    fig = plt.subplot(2,2,1);
    M = convertGraphDicToSparseMatrix(gVpSrp)
    plt.spy(M, aspect='auto',marker='s',markersize=5)
    plt.title("gVpSrp")
    fig = plt.subplot(2,2,2);
    M = convertGraphDicToSparseMatrix(gVpStp)
    plt.spy(M, aspect='auto',marker='s',markersize=5)
    plt.title("gVpStp")
    fig = plt.subplot(2,2,3);
    M = convertGraphDicToSparseMatrix(gSrp)
    plt.spy(M, aspect='auto',marker='s',markersize=5)
    plt.title("gSrp")
    fig = plt.subplot(2,2,4);
    M = convertGraphDicToSparseMatrix(gStp)
    plt.spy(M, aspect='auto',marker='s',markersize=5)
    plt.title("gStp")

    # plt.xticks(np.arange(gSpPrint[numResidualPtsSp-1][0]+1))
    # plt.yticks(np.arange(gVpPrint[numResidualPtsVp-1][0]+1))

  # -----------------------------------------------------
  # number of pts where we store the state is the same as residual points
  numStatePtsVp = numResidualPtsVp
  numStatePtsSrp = numResidualPtsSrp
  numStatePtsStp = numResidualPtsStp
  numStatePts   = numStatePtsVp + numStatePtsSrp + numStatePtsStp
  print ("numStatePtsVp = ", numStatePtsVp)
  print ("numStatePtsSrp = ", numStatePtsSrp)
  print ("numStatePtsStp = ", numStatePtsStp)
  print ("numStatePts = ", numStatePts,
         " which is = ", numStatePts/totNumDofs*100, " % of full mesh")


  printMeshInfoFileThreeVarGroups("full", nth, nr,
                                  thLDeg, thRDeg, rCmb, rSurf, dth, dr,
                                  numStatePtsVp, numResidualPtsVp,
                                  numStatePtsSrp, numResidualPtsSrp,
                                  numStatePtsStp, numResidualPtsStp)

  writeConnectivityFileFullMesh(gVpSrp, coordsVp, "graph_vpsrp.dat")
  writeConnectivityFileFullMesh(gVpStp, coordsVp, "graph_vpstp.dat")
  writeVpStencilCoeffsFileFullMesh(coeffVpSrp, "coeff_vpsrp.dat")
  writeVpStencilCoeffsFileFullMesh(coeffVpStp, "coeff_vpstp.dat")
  writeConnectivityFileFullMesh(gSrp, coordsSrp, "graph_srp.dat")
  writeConnectivityFileFullMesh(gStp, coordsStp, "graph_stp.dat")

  # move all files to target directory
  os.system("mv *.dat {}".format(workDir))

  # -----------------------------------------------------
  if plotting != "none":
    fig = plt.figure(0)
    ax = fig.add_subplot(111, projection='polar')
    ax.set_rlabel_position(180)
    ax.set_aspect(1.0)
    plotEarthSurf(ax)
    plotCMB(ax)

    th_fm = np.asarray([v[0] for (k,v) in coordsVp.items()])
    r_fm  = np.asarray([v[1] for (k,v) in coordsVp.items()])
    gids  = np.asarray([k for (k,v) in coordsVp.items()])
    ax.scatter(th_fm, r_fm, marker='o', s=4, edgecolor='r', facecolor='r')
    for i in range(len(gids)):
      ax.text(th_fm[i], r_fm[i], str( gids[i] ),
               verticalalignment='top', horizontalalignment='center', fontsize=10, color='r')

    th_fm = np.asarray([v[0] for (k,v) in coordsSrp.items()])
    r_fm  = np.asarray([v[1] for (k,v) in coordsSrp.items()])
    gids  = np.asarray([k for (k,v) in coordsSrp.items()])
    ax.scatter(th_fm, r_fm, marker='s', s=4, edgecolor='k', facecolor='k')
    for i in range(len(gids)):
      ax.text(th_fm[i], r_fm[i], str( gids[i] ),
               verticalalignment='top', horizontalalignment='center', fontsize=10, color='k')

    th_fm = np.asarray([v[0] for (k,v) in coordsStp.items()])
    r_fm  = np.asarray([v[1] for (k,v) in coordsStp.items()])
    gids  = np.asarray([k for (k,v) in coordsStp.items()])
    ax.scatter(th_fm, r_fm, marker='^', s=4, edgecolor='b', facecolor='b')
    for i in range(len(gids)):
      ax.text(th_fm[i], r_fm[i], str( gids[i] ),
               verticalalignment='top', horizontalalignment='center', fontsize=10, color='b')
    plt.show()
