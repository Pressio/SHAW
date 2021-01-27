
import matplotlib.pyplot as plt
import sys, os, time
import numpy as np
from numpy import linspace, meshgrid
import collections

class UnitCell:
  # the unit building block for the mesh
  def __init__(self, cellGID, cell_ij, numDof, labels,
               cell_thCoords, cell_rrCoords, onSymmetryAxis):
    # the cellGID is the global ID identifier of the cell within the grid
    self.cellGID_ = cellGID
    # the i,j indices of this cell (this is wrt to all the cells in domain)
    self.ij_ = cell_ij
    # num of dofs inside this cell
    self.numDof_ = numDof
    # which dofs are defined inside this cell
    self.dofLabels_ = labels
    # the theta coords of all grid points inside this cell
    self.th_ = cell_thCoords
    # the r coords of all grid points inside this cell
    self.rr_ = cell_rrCoords
    self.markers = {'vp': 'o', 'stp': 's', 'srp': 's'}
    self.edgeCol = {'vp': 'r', 'stp': 'k', 'srp': 'k'}
    self.faceCol = {'vp': 'r', 'stp': 'k', 'srp': 'none'}
    # the GIDs of the dofs I own
    self.myDofGIDs_ = []
    # if this is a cell next to symmetry axes
    self.onSymmetryAxis_ = onSymmetryAxis

  def isOnSymmetryAxis(self):
    return self.onSymmetryAxis_

  def getCellGID(self):
    return self.cellGID_

  def getDofIndex(self, targetLabel):
    assert( targetLabel == 'vp' or targetLabel=='stp' or targetLabel=='srp')
    return self.dofLabels_.index(targetLabel)

  def dofExists(self, targetLabel):
    return targetLabel in self.dofLabels_

  def getDofGID(self, targetLabel):
    assert( targetLabel == 'vp' or targetLabel=='stp' or targetLabel=='srp')
    index = self.dofLabels_.index(targetLabel)
    return self.myDofGIDs_[index]

  def getDofCoords(self, targetLabel):
    assert( targetLabel == 'vp' or targetLabel=='stp' or targetLabel=='srp')
    index = self.dofLabels_.index(targetLabel)
    return [self.th_[index], self.rr_[index]]

  def getCellIJ(self):
    return self.ij_

  def computeDofsGIDsStartingFrom(self, startGIDvp, startGIDsp):
    for i in range(self.numDof_):
      if self.dofLabels_[i] == 'vp':
        startGIDvp+=1
        self.myDofGIDs_.append(startGIDvp)
      else:
        startGIDsp+=1
        self.myDofGIDs_.append(startGIDsp)
    return [startGIDvp, startGIDsp]

  def computeDofsGIDsSplitVarsStartingFrom(self, startGIDvp, startGIDsrp, startGIDstp):
    for i in range(self.numDof_):
      if self.dofLabels_[i] == 'vp':
        startGIDvp+=1
        self.myDofGIDs_.append(startGIDvp)
      elif self.dofLabels_[i] == 'srp':
        startGIDsrp+=1
        self.myDofGIDs_.append(startGIDsrp)
      else:
        startGIDstp+=1
        self.myDofGIDs_.append(startGIDstp)
    return [startGIDvp, startGIDsrp, startGIDstp]

  def plotPoints(self, ax):
    for i in range(self.numDof_):
      label = self.dofLabels_[i]
      ax.scatter(self.th_[i], self.rr_[i], marker=self.markers[label], s=2,
                 facecolors=self.faceCol[label], edgecolors=self.edgeCol[label])
      #ax.text(self.th_[i], self.rr_[i], str(self.myDofGIDs_[i]),
      #        verticalalignment='top', horizontalalignment='center', fontsize=14,
      #        color=self.edgeCol[label])
