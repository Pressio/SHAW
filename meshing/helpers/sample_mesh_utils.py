#!/usr/bin/env python

import sys, os, time
import numpy as np
from numpy import linspace, meshgrid
import collections

#------------------------------------------------------------------------------------------
def buildGraphSampleMesh(fullMeshCells, nth, nr, targetPct, splitVars = True):
  totalNumOfCells = len(fullMeshCells)

  # create a list of cells IDs we want for sample mesh
  targetListOfCells = createRandomListOfGIDsForSampleMesh(totalNumOfCells, targetPct)

  # indicesVp = np.loadtxt("sm_indices_vp.dat")
  # indicesSp = np.loadtxt("sm_indices_sp.dat")
  # #targetListOfCells = []
  # for k in fullMeshCells:
  #   cellID = k.getCellGID()
  #   vpGID = k.getDofGID('vp')
  #   if vpGID in indicesVp:
  #     targetListOfCells.append(cellID)

  #   if k.dofExists('srp'):
  #     srpGID = k.getDofGID('srp')
  #     if srpGID in indicesSp:
  #       targetListOfCells.append(cellID)

  #   if k.dofExists('stp'):
  #     stpGID = k.getDofGID('stp')
  #     if stpGID in indicesSp:
  #       targetListOfCells.append(cellID)

  # targetListOfCells = np.array([0]) #np.arange(0, totalNumOfCells, 20)
  for k in fullMeshCells:
    cellID = k.getCellGID()
    vpGID = k.getDofGID('vp')
    if vpGID==27665:
      targetListOfCells = np.append(targetListOfCells, cellID)
    vpTheta = np.degrees( k.getDofCoords('vp')[0] )
    # if vpTheta > 24 and vpTheta < 26:
    #   targetListOfCells = np.append(targetListOfCells, cellID)
    # vprr = k.getDofCoords('vp')[1]
    # if vprr > 6371.0-710 and vprr < 6371.-700:
    #   targetListOfCells = np.append(targetListOfCells, cellID)

  targetListOfCells = np.unique(targetListOfCells)
  # sort for convenience
  targetListOfCells.sort()

  # call impl
  return buildGraphImpl(fullMeshCells, nth, nr, targetListOfCells, splitVars)
