#!/usr/bin/env python

import matplotlib.pyplot as plt
import sys, os, time
import numpy as np
from numpy import linspace, meshgrid
from matplotlib import cm
import collections
from argparse import ArgumentParser
import random
import scipy.sparse as sp
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import reverse_cuthill_mckee

sys.path.insert(0, './helpers')
#thisDir = os.getcwd()
#sys.path.append(os.path.abspath(thisDir+"/helpers"))

from earth_params import earthRadius, mantleThickness, cmbRadius
from fullMeshWithThreeVarGroups import *
from fullMeshWithTwoVarGroups import *

#--------------------------------------------------------------
def main(workDir, nth, nr, samplingType, targetPct, plotting,
         thL, thLDeg, thR, thRDeg, rSurf, rCmb, debugPrint):

  # dimensions along r and theta
  L = [thR-thL, rSurf-rCmb]
  print ("Bounds for th = ", thL, thR,    " [rad] with L= ", L[0], " [rad]")
  print ("Bounds for r  = ", rCmb, rSurf, " [km] with  L= ", L[1], " [km]")

  # compute cell size
  dth, dr = L[0]/(nth-1), L[1]/(nr-1)
  print ("Cell spacing: dth = ", dth, " [rad] ")
  print ("Cell spacing: dr = ", dr, " [km]")

  print ("Angular spacing at CMB [km]: ", dth * cmbRadius)
  print ("Angular spacing at earth surface [km]: ", dth * earthRadius)

  #-----------------
  #--- FULL MESH ---
  #-----------------
  varGroups=2
  if (samplingType=="full"):
    if (varGroups==2):
      mainFullMeshWithTwoVarGroups(workDir, nth, nr, plotting,
                                   thL, thLDeg, thR, thRDeg,
                                   rSurf, rCmb, debugPrint,
                                   L, dth, dr)
    elif (varGroups==3):
      print("Fully splitting the dofs three-ways not supported yet")
      sys.exit(1)
      # #mainFullMeshWithThreeVarGroups(workDir, nth, nr, plotting,
      #                                thL, thLDeg, thR, thRDeg,
      #                                rSurf, rCmb, debugPrint,
      #                                L, dth, dr)

  #-------------------
  #--- Sample MESH ---
  #-------------------
  elif (samplingType=="random"):
    print("Sample mesh implementation not yet deloyed")
    sys.exit(1)

  #-----------------
  #--- not valid ---
  #-----------------
  else:
    print("unknown value for sample-type")
    sys.exit(1)


###############################
if __name__== "__main__":
###############################
  parser = ArgumentParser()
  parser.add_argument("-nr", "--nr",    type=np.int64, dest="nr")
  parser.add_argument("-nth", "--nth",  type=np.int64, dest="nth")

  parser.add_argument("-working-dir", "--working-dir",
                      dest="workDir", default=".",
                      help="Parent dir where I will create a subdir with all generated mesh files to.")

  parser.add_argument("-thLeft", "--thLeft",
                      type=float, dest="thL",
                      default=0.,
                      help="The theta coord in deg of the left bound along theta, default = 0")
  parser.add_argument("-thRight", "--thRight",
                      type=float, dest="thR",
                      default=180.,
                      help="The theta coord in deg of the right bound along theta, default = 180")

  parser.add_argument("-rCmb", "--rCmb",
                      type=float, dest="rCmb",
                      default=cmbRadius,
                      help="The r coordinate [km] of the CMB, default = " + str(cmbRadius))

  parser.add_argument("-rSurf", "--rSurf",
                      type=float, dest="rSurf",
                      default=earthRadius,
                      help="The r coordinate [km] of the earth surface, default = " + str(earthRadius))

  parser.add_argument("-plotting", "--plotting",
                      dest="plotting",
                      default="none",
                      help="What type of plotting you want:\n"+
                           "use <show> for showing plots,\n"+
                           "use <print> for printing only, \n"+
                           "use <none> for no plots")

  parser.add_argument("-sampling-type", "--sampling-type",
                      dest="samplingType",
                      default="full",
                      help="What type of mesh you need:\n"+
                           "use <full> for creating connectivity for full mesh,\n"+
                           "use <random> for creating connectivity for sample mesh")

  # parser.add_argument("-var-groups", "--var-groups",
  #                     dest="varGroups", type=int,
  #                     default=2,
  #                     help="How to split and group variables:\n"+
  #                          "use <2> for two groups: vp and sp=(srp,stp),\n"+
  #                          "use <3> for fully splitting vp,srp,stp")

  parser.add_argument("-debug-print", "--debug-print",
                      type=int, dest="debugPrint",
                      default=0,
                      help="If you want debug prints, disabled by default\n")

  parser.add_argument("-pct", "--pct",
                      dest="targetPct",
                      default=-1,
                      help="The target % of elements where to compute residual.\n"+
                      "It must be 0 < pct <= 100.\n"+
                      "Note that more points are needed for the state\n"+
                      "because for every residual location I need to account for neighbors."+
                      "You do not need this if you select -sampling-type=full,\n"+
                      "since in that case we sample the full mesh.\n"+
                      "But you need to set it for -sampling-type=random")

  ##################
  # parse all args
  ##################
  args = parser.parse_args()

  # check for valid inputs of grid sizes
  if (args.nr <= 0): print(" --nr must be > 0"); sys.exit(1)
  if (args.nth <= 0): print(" --nth must be > 0"); sys.exit(1)

  # check that domain limits make sense
  if (args.rCmb >= args.rSurf):
    print(" rCmb must be < rSurf, i.e. radius of CMB must be smaller than the radius of earth surface")
    sys.exit(1)
  if (args.rCmb <= 0):
    print(" rCmb must be > 0, i.e. radius of CMB must be greater than zero")
    sys.exit(1)
  if (args.rSurf <= 0):
    print(" rSurf must be > 0")
    sys.exit(1)

  if (args.thL >= args.thR):
    print(" thLeft must be < thRight ")
    sys.exit(1)

  # convert from degress to rad the theta limits
  thLDeg, thRDeg = args.thL, args.thR
  args.thL = args.thL * np.pi / 180.
  args.thR = args.thR * np.pi / 180.

  if (args.samplingType == "random"):
    if (float(args.targetPct) <= 0 or float(args.targetPct) > 100):
      print(" --pct must be 0 < --pct <= 100")
      sys.exit(1)

  # fix seed for reproducibility
  random.seed( 5243543 )

  main(args.workDir, args.nth, args.nr, args.samplingType,
       float(args.targetPct), args.plotting,
       float(args.thL), float(thLDeg),
       float(args.thR), float(thRDeg),
       float(args.rSurf), float(args.rCmb),
       args.debugPrint)
