#!/usr/bin/env python

import sys, os, time
import numpy as np
import collections
from unit_cell import UnitCell

def createCellsFullMesh(nth, nr, thL, thR, rSurf, rCmb, dth, dr, varGroups = 3):
  # origin of the full grid
  gridO = [thL, rCmb]
  # list of all cells
  myCells = []

  # total # of velocity points
  totVelPts = nr*nth

  # loop and create each cell by moving along theta first, then to next radius.
  # This leads to a natural orderning of the cells
  # j: indexing the cells along r
  # i: indexing the cells along theta

  # cellID enumerates each cell wrt the whole domain
  cellID = 0

  # set counters to track the global ID of every dof in the system
  # (need to start from -1 so that first dof has GID=0)
  vp_dofGID = -1
  if varGroups == 3:
    srp_dofGID = -1
    stp_dofGID = -1
  elif varGroups == 2:
    sp_dofGID = -1

  for j in range(nr):
    for i in range(nth):
      # ever cell has at least 1 dof, which is the velocity vp
      thisCellNumDofs = 1

      # all cells have velocity dof
      thisCellLabels = ['vp']

      # this is the index pair identifying the cell wrt the 2d grid
      cell_ij = [i,j]

      # find origin of current cell: it corresponds to the velocity grid point
      # Express it as a shift wrt the origin of the full domain
      cellOrigin = [gridO[0] + i*dth, gridO[1] + j*dr]

      # all cells have at least theta,r for vp point which is also origin of the cell
      vpTh, vpR = cellOrigin[0], cellOrigin[1]
      thisCellth = [vpTh]
      thisCellrr = [vpR]

      # all cells below the earth surface also have the srp dof
      if j < nr-1:
        # increment the num of dofs inside this cell
        thisCellNumDofs += 1
        # append the label for srp
        thisCellLabels.append('srp')
        # srp point has same theta as vp, but r is shifted upwards by dr/2
        thisCellth.append(vpTh)
        thisCellrr.append(vpR+dr*0.5)

      # if i < nth-1, the cells also have the stp dof
      if i < nth-1:
        thisCellNumDofs += 1
        thisCellLabels.append('stp')
        # stp point has same r as vp, but theta is shifted by dth/2
        thisCellth.append(vpTh + dth*0.5)
        thisCellrr.append(vpR)

      onSymmetryAxis = True if (i==0 or i==nth-1) else False
      # create and append new cell
      myCells.append(UnitCell(cellID, cell_ij,
                              thisCellNumDofs, thisCellLabels,
                              thisCellth, thisCellrr, onSymmetryAxis))

      # pass to cell object the current global ID enumerating all the dofs,
      # ask the cell object to compute the GIDs inside and return the updated count
      if varGroups == 3:
        [vp_dofGID, srp_dofGID, stp_dofGID] = myCells[-1].computeDofsGIDsSplitVarsStartingFrom(vp_dofGID,
                                                                                               srp_dofGID,
                                                                                               stp_dofGID)
      elif varGroups == 2:
        [vp_dofGID, sp_dofGID] = myCells[-1].computeDofsGIDsStartingFrom(vp_dofGID, sp_dofGID)

      # increment the cell indexing
      cellID+=1
  return myCells


#------------------------------------------------------------------------------------------
def storeCoordsAndLabelsFullMesh(fullMeshCells, nth, nr, varGroups = 3):
  # the total number of cells we have
  numCells = len(fullMeshCells)

  # to store coordinates, use dic:
  # key = the global ID of a dof
  # value = a list with [theta, rad] of the odf
  coordsVp = {}
  isOnSymVp = {}
  if varGroups == 2:
    coordsSp = {}
  else:
    coordsSrp = {}
    coordsStp = {}

  if varGroups == 2:
    # dofLabels is a map such that:
    # key = the global ID of a dof
    # value = the label of the odf, we use 1=srp, 2=stp
    # this is only needed when we DO not split srp,stp
    dofLabelsSp = {}
    # to flag points on symmetry axes
    isOnSymSp = {}

  # --- loop over list of target cell IDs we want ---
  for k in range(numCells):
    iCell = fullMeshCells[k]
    # find the ij of this cell
    cellI, cellJ = iCell.getCellIJ()[0], iCell.getCellIJ()[1]

    #----------------------
    # --- deal with vp ---
    #----------------------
    # get GID of the vp dof at current cell
    vpGID = iCell.getDofGID('vp')
    coordsVp[vpGID] = iCell.getDofCoords('vp')
    isOnSymVp[vpGID] = iCell.isOnSymmetryAxis()

    #----------------------
    # --- deal with srp ---
    #----------------------
    # srp exists only for cells with j<=nr-2
    if cellJ <= nr-2:
      label = 'srp'
      srpGID = iCell.getDofGID(label)
      if varGroups==3:
        coordsSrp[srpGID] = iCell.getDofCoords(label)
      else:
        coordsSp[srpGID] = iCell.getDofCoords(label)
        dofLabelsSp[srpGID] = 1
        isOnSymSp[srpGID] = iCell.isOnSymmetryAxis()

    #----------------------
    # --- deal with stp ---
    #----------------------
    # stp exists only for cells with i<=nth-2
    if cellI <= nth-2:
      label = 'stp'
      stpGID = iCell.getDofGID(label)
      if varGroups==3:
        coordsStp[stpGID] = iCell.getDofCoords(label)
      else:
        coordsSp[stpGID] = iCell.getDofCoords(label)
        dofLabelsSp[stpGID] = 2
        # remember that only srp points can be on symmetry axis, all stp are NOT
        isOnSymSp[stpGID] = 0

  if varGroups==3:
    return [coordsVp, coordsSrp, coordsStp]
  else:
    return [coordsVp, coordsSp, dofLabelsSp, isOnSymVp, isOnSymSp]


# #------------------------------------------------------------------------------------------
# def buildGraphFullMesh(fullMeshCells, nth, nr):
#   # here we want to use all cells, so we call the buildGraphImpl
#   # passing a list of all cell IDs since we consider all of them, not just a subset
#   # max ID of the cells is equal to nr*nth-1
#   targetListOfCells = range(0, nth*nr)
#   # call impl
#   return buildGraphImpl(fullMeshCells, nth, nr, targetListOfCells, False)

# #------------------------------------------------------------------------------------------
# def buildGraphFullMeshSplitVars(fullMeshCells, nth, nr):
#   # here we want to use all cells, so we call the buildGraphImpl
#   # passing a list of all cell IDs since we consider all of them, not just a subset
#   # max ID of the cells is equal to nr*nth-1
#   targetListOfCells = range(0, nth*nr)
#   # call impl
#   return buildGraphImpl(fullMeshCells, nth, nr, targetListOfCells, True)
