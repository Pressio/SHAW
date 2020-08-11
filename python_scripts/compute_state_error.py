#!/usr/bin/env python

import numpy as np
import sys, re, os
import matplotlib.pyplot as plt
from argparse import ArgumentParser
from numpy import linalg as la

#=========================================
def str2bool(v):
  if isinstance(v, bool):
    return v
  if v.lower() in ('yes', 'true', 't', 'y', '1'):
    return True
  elif v.lower() in ('no', 'false', 'f', 'n', '0'):
    return False
  else:
    raise argparse.ArgumentTypeError('Boolean value expected.')

#=========================================
def loadStates(fomFile, approxFile, dofName):
  # load states
  if not os.path.exists(fomFile):
    print("fom final state {} does not exist".format(fomFile))
    sys.exit(1)

  if not os.path.exists(approxFile):
    print("approx final state {} does not exist".format(approxFile))
    sys.exit(1)

  print("reading fom    state {}".format(fomFile))
  print("reading approx state {}".format(approxFile))

  # load data (skip first row because it contains the size)
  fomState              = np.loadtxt(fomFile, skiprows=1)
  fomStateReconstructed = np.loadtxt(approxFile, skiprows=1)
  fomNRows = fomState.shape[0]
  approxNRows = fomStateReconstructed.shape[0]
  # if things are correct, this should be a single vector and approx/fom match sizes
  assert( fomNRows == approxNRows )

  return [fomState, fomStateReconstructed]

#=========================================
def computeErrors(fomState, fomStateReconstructed, dofName):
  error = fomState-fomStateReconstructed
  print(" {}_fom_minmax    = {} {}".format(dofName, np.min(fomState), np.max(fomState)))
  print(" {}_approx_minmax = {} {}".format(dofName, np.min(fomStateReconstructed), np.max(fomStateReconstructed)))
  print(" {}_err_minmax    = {} {}".format(dofName, np.min(error), np.max(error)))

  fomL2Norm, fomLinfNorm = la.norm(fomState), la.norm(fomState, np.inf)
  errL2Norm   = [ la.norm(error),         la.norm(error)/fomL2Norm ]
  errLinfNorm = [ la.norm(error, np.inf), la.norm(error, np.inf)/fomLinfNorm ]
  print(" {}_err_abs_rel_ltwo_norms = {} {}".format(dofName, errL2Norm[0],   errL2Norm[1]))
  print(" {}_err_abs_rel_linf_norms = {} {}".format(dofName, errLinfNorm[0], errLinfNorm[1]))


###############################
if __name__== "__main__":
###############################
  parser = ArgumentParser()
  parser.add_argument("-dryrun", "--dryrun", "-dr", "--dr",
                      dest="dryRun", type=str2bool, default=True,
                      help="True: creates directory structures/files, does not run. Default=True.")

  parser.add_argument("-fom-state", "--fom-state", "-fomstate", "--fomstate",
                      dest="fomState", default="empty",
                      help="Full path to fom state. Must be set.")

  parser.add_argument("-approx-state", "--approx-state", "-approxstate", "--approxstate",
                      dest="approxState", default="empty",
                      help="Full path to approx state. Must be set.")

  parser.add_argument("-dof-name", "--dof-name", "-dofname", "--dofname",
                      dest="dofName", default="empty",
                      help="Which dof: vp/sp to compute error for.")

  # parse args
  args = parser.parse_args()
  assert(args.fomState != "empty")
  assert(args.approxState != "empty")
  assert(args.dofName != "empty")

  dofName = args.dofName
  [fomState, fomStateReconstructed] = loadStates(args.fomState, args.approxState, dofName)
  computeErrors(fomState, fomStateReconstructed, dofName)








  # error = fomState - fomStateReconstructed
  # fomBounds = [np.min(fomState), np.max(fomState)]
  # nr, nth = 256, 1024
  # cc = np.loadtxt(args.fomDir+"/coords_" + dofName + ".txt")
  # th, r = cc[:,0], cc[:, 1]
  # th, r = th.reshape((nr,nth)), r.reshape((nr,nth))

  # fig1 = plt.figure(1)
  # ax1 = fig1.add_subplot(111, projection='polar')
  # h1=ax1.pcolormesh(th, r, fomState.reshape((nr,nth)), cmap="twilight_shifted", shading = "flat",
  #                   vmin=fomBounds[0], vmax=fomBounds[1])
  # ax1.set_rlabel_position(260)
  # fig1.colorbar(h1)

  # fig2 = plt.figure(2)
  # ax2 = fig2.add_subplot(111, projection='polar')
  # h2=ax2.pcolormesh(th, r, error.reshape((nr,nth)), cmap="binary", shading = "flat",
  #                   vmin=fomBounds[0], vmax=fomBounds[1]) #vmin=np.min(error), vmax=np.max(error))
  # ax2.set_rlabel_position(260)
  # fig2.colorbar(h2)

  # fig3 = plt.figure(3)
  # ax3 = fig3.add_subplot(111, projection='polar')
  # h3=ax3.pcolormesh(th, r, fomStateReconstructed.reshape((nr,nth)), cmap="twilight_shifted", shading = "flat",
  #   vmin=fomBounds[0], vmax=fomBounds[1])
  # ax3.set_rlabel_position(260)
  # fig3.colorbar(h3)

  # plt.show()
