
import sys, os, time
import numpy as np
from numpy import linspace, meshgrid
import random

def printDicPretty(d):
  for key, value in d.items():
    print(str(key), value)

def printDicPretty2(d1, d2):
  for (k1, v1), (k2, v2) in zip(d1.items(), d2.items()):
    print(str(k1), v1, " -- ", str(k2), v2)


def createRandomListOfGIDsForSampleMesh(numCells, targetPct):
  # assert that the percentage is given in specific format, not as a fraction
  assert( targetPct > 0 and targetPct<=100 )

  # convert targetPct (which is between 0 and 100) to fraction
  targetPctFrac = targetPct * 1e-2

  # number of cells we want in the sample mesh
  targetSMSize = targetPctFrac * numCells

  # need to get the floor to get an integer
  targetSMSize = np.int64(np.floor(targetSMSize))

  # get unique list of random IDs picked from the full IDs
  rGIDs = random.sample(range(numCells), targetSMSize)

  return rGIDs
