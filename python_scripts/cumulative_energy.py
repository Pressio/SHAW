#!/usr/bin/env python

import numpy as np
import sys, re
from argparse import ArgumentParser

def computeCumE(s, target):
  sSq = np.square(s)
  den = np.sum(sSq)
  rsum = 0.
  for i in range(0, len(s)):
    rsum += sSq[i]
    ratio = (rsum/den)
    if ratio >= target:
      return i
  return len(s)

if __name__== "__main__":
  parser = ArgumentParser()
  parser.add_argument("-working-dir", "--working-dir", "-wdir", "--wdir",
                      dest="workDir", default="empty",
                      help="Target dir where to work. Must be set.")

  parser.add_argument("-percent", "--percent",
                      dest="pct", default=999, type=np.float,
                      help="Target fraction to keep")

  # parse all args
  args = parser.parse_args()
  assert(args.workDir != "empty")

  # convert percentage to decimal
  target = float(args.pct)/100.
  print("target {}".format(target))

  # load data
  sVp = np.loadtxt(args.workDir+"/sva_vp")
  sSp = np.loadtxt(args.workDir+"/sva_sp")

  # compute cumulative energy
  if (target == 1.):
    nVp, nSp = len(sVp), len(sSp)
    print ("Nvp: {} Nsp: {}".format(nVp, nSp) )
  else:
    nVp = computeCumE(sVp, target)
    nSp = computeCumE(sSp, target)

  print ("Nvp: {} Nsp: {}".format(nVp, nSp) )
  maxN = max(nVp, nSp)
  print ("numBasis: {}".format(maxN) )
