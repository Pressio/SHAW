#!/usr/bin/env python

import sys, os, time
import numpy as np
import collections

def buildGraphImpl(fullMeshCells, nth, nr, targetListOfCells, varGroups):
  # the total number of cells we have
  numCells = len(fullMeshCells)

  # graph is an adjecency list (here we use a dictionary) for each dof:
  # - key   = the global ID of a dof
  # - value = list of global IDs for the dofs needed at west, north, eest, south
  if varGroups==2:
    graphVp = collections.OrderedDict()
    graphSp = collections.OrderedDict()
    # coeffVp contains the sign (-1 or 1) to identify where a Vp is located
    # such that we can handle the boundaries easily
    # - key   = the global ID of a dof
    # - value = value needed for the west, north, eest, south dof
    coeffVp = collections.OrderedDict()

  else:
    graphVpSrp = collections.OrderedDict()
    graphVpStp = collections.OrderedDict()
    coeffVpSrp = collections.OrderedDict()
    coeffVpStp = collections.OrderedDict()
    graphSrp = collections.OrderedDict()
    graphStp = collections.OrderedDict()

  # --- loop over list of target cell IDs we want ---
  for k in targetListOfCells:
    iCell = fullMeshCells[k]

    # find the ij of this cell
    cellI, cellJ = iCell.getCellIJ()[0], iCell.getCellIJ()[1]

    # temporary list of neighbors for a vp points listed as west, north, east, south (order matters)
    tmpListVp  = np.zeros(4, dtype=np.int64)
    # tmp fro storing multiplicative coeff for vp stencil
    tmpCoeffVp = np.zeros(4, dtype=float)
    # stp points only have west, east (order matters)
    tmpListStp = np.zeros(2, dtype=np.int64)
    # srp points only have north, south (order matters)
    tmpListSrp = np.zeros(2, dtype=np.int64)

    #****************************************
    #           deal with vp
    #
    # vp exists everywhere but
    # need to be careful near the boundaries
    #****************************************

    # get GID of the vp dof at current cell
    vpGID = iCell.getDofGID('vp')

    #-----------------------------------
    # if cell is NOT close to boundary at theta=thetaLeft, stp dof on the left exists
    if cellI >=1:
      # get the gid of west stp, which is in the cell on the left
      # since we use natural order for cell enumerating, the cell on the left is at k-1
      westCell = fullMeshCells[k-1]
      tmpListVp[0]  = westCell.getDofGID('stp')
      tmpCoeffVp[0] = 1.
    # if cell is at theta=thetaLeft
    if cellI==0:
      # here we are at thetaLeft, and we know we need the stp on the left which
      # does not exits. but we also know that the stp on the left is set
      # to be equal to negative the value of
      # stp on the right so we use trick of setting the west GID to be the negative of
      # the stp dof on the right
      tmpListVp[0]  = iCell.getDofGID('stp')
      tmpCoeffVp[0] = 0.

    #-----------------------------------
    # if cell is NOT close to thetaRight
    if cellI<nth-1:
      # get the gid of east stp
      tmpListVp[2]  = iCell.getDofGID('stp')
      tmpCoeffVp[2] = 1.
    if cellI==nth-1:
      # here we are right at thetaRight, and we know we need the stp on the right
      # which does not exits. but we also know that the stp on the right is set to
      # be equal to negative the stp on the left so we set the east GID to be
      # the negative of the stp dof on the west
      westCell = fullMeshCells[k-1]
      tmpListVp[2]  = westCell.getDofGID('stp')
      tmpCoeffVp[2] = 0.

    #-----------------------------------
    # if cell is NOT at earth surface
    if cellJ<nr-1:
      # get the gid of the north srp
      tmpListVp[1]  = iCell.getDofGID('srp')
      tmpCoeffVp[1] = 1.
    # if cell is at earth surface
    if cellJ==nr-1:
      northCell = fullMeshCells[k-nth]
      tmpListVp[1]  = northCell.getDofGID('srp')
      tmpCoeffVp[1] = -1.

    #-----------------------------------
    # if cell is NOT at CMB
    if cellJ>=1:
      # get the gid of south srp
      southwCell = fullMeshCells[k-nth]
      tmpListVp[3]  = southwCell.getDofGID('srp')
      tmpCoeffVp[3] = 1.
    # if cell is right next to cmb
    if cellJ==0:
      tmpListVp[3]  = iCell.getDofGID('srp')
      tmpCoeffVp[3] = -1.


    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # at the symmetry axes, set all coeffs = 0
    # since v undefined here, this also means that
    # srp will be zero here because srp uses the v on this axess
    if cellI==0 or cellI==nth-1:
      tmpCoeffVp = np.zeros(4)
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    # add to graph
    if (varGroups==2):
      graphVp[vpGID] = tmpListVp
      coeffVp[vpGID] = tmpCoeffVp
    else:
      # for vpSrp, only take north and south elements
      graphVpSrp[vpGID] = tmpListVp[1::2]
      coeffVpSrp[vpGID] = tmpCoeffVp[1::2]
      # for vpStp, only take west and east elements
      graphVpStp[vpGID] = tmpListVp[0::2]
      coeffVpStp[vpGID] = tmpCoeffVp[0::2]

    #***********************************************
    #           --- deal with srp ---
    #***********************************************
    # srp exists only for cells with j<=nr-2
    if cellJ <= nr-2:
      srpGID = iCell.getDofGID('srp')
      # get the cell that is north
      northCell = fullMeshCells[k+nth]
      # needs the north velocity
      tmpListSrp[0] = northCell.getDofGID('vp')
      # needs south velocity
      tmpListSrp[1] = iCell.getDofGID('vp')
      if varGroups==2:
        graphSp[srpGID] = tmpListSrp
      else:
        graphSrp[srpGID] = tmpListSrp

    #***********************************************
    #           --- deal with stp ---
    #***********************************************
    # stp exists only for cells with i<=nth-2
    if cellI <= nth-2:
      stpGID = iCell.getDofGID('stp')
      # needs west velocity
      tmpListStp[0] = iCell.getDofGID('vp')
      # get the cell that is neighbor on right
      eastCell = fullMeshCells[k+1]
      # needs east velocity
      tmpListStp[1] = eastCell.getDofGID('vp')
      if varGroups==2:
        graphSp[stpGID] = tmpListStp
      else:
        graphStp[stpGID] = tmpListStp

  if varGroups==2:
    return [graphVp, coeffVp, graphSp]
  else:
    return [graphVpSrp, graphVpStp, coeffVpSrp, coeffVpStp, graphSrp, graphStp]
