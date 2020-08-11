#!/usr/bin/env python

import numpy as np
import sys, re, os, yaml
import matplotlib.pyplot as plt
from argparse import ArgumentParser
from numpy import linalg as la
from scipy.interpolate import interp1d, CubicSpline
from scipy.interpolate import UnivariateSpline, Rbf

from utils import str2bool

#=========================================
def loadStates(fomDir, romDir, dofName):
  fomFile = fomDir+"/finalFomState_"+dofName+"_0"
  romFile = romDir+"/finalFomState_"+dofName+"_0"

  # load fom and rom final state
  if not os.path.exists(fomFile):
    print("fom final state {} does not exist".format(fomFile))
    sys.exit(1)

  if not os.path.exists(romFile):
    print("rom final state {} does not exist".format(romFile))
    sys.exit(1)

  # load data (skip first row because it contains the size)
  fomState              = np.loadtxt(fomFile, skiprows=1)
  fomStateReconstructed = np.loadtxt(romFile, skiprows=1)
  fomNRows = fomState.shape[0]
  romNRows = fomStateReconstructed.shape[0]
  # if things are correct, this should be a single vector and rom/fom match sizes
  assert( fomNRows == romNRows )
  return [fomState, fomStateReconstructed]

#=========================================
def getListOfDirs(dryRun, scenario, workDir, key, sortByRunID):
  # find all dirs with "train" in the name to use as trains
  if key == "train":
    dirsFullPath = [workDir+'/'+d for d in os.listdir(workDir) if key in d]
  elif key=="test":
    dirsFullPath = [workDir+'/'+d for d in os.listdir(workDir) if key in d and "fom" in d]
  else:
    sys.exit("Invalid key value for getListOfDirs")

  if sortByRunID:
    # sort based on the ID, so need to extract ID which is last of dir name
    def func(elem): return int(elem.split('_')[-1])
    return sorted(dirsFullPath,key=func)
  else:
    dirsFullPath

#=========================================
def extractTargetValue(dr, dirPath, scenario):
  with open(dirPath+'/input.yaml') as file:
    # The FullLoader parameter handles the conversion from YAML
    # scalar values to Python the dictionary format
    ifile = yaml.load(file, Loader=yaml.FullLoader)
    if scenario == 1:
      return float( ifile['material']['layer2']['velocity'][0] )
    elif scenario == 2:
      return float( ifile['source']['signal']['period'] )
    else:
      sys.exit("Invalid scenario")

#=========================================
def getTargetValues(dr, scenario, dirs):
  values = np.zeros(len(dirs))
  for i,iDir in enumerate(dirs):
    values[i] = extractTargetValue(dr, iDir, scenario)
  return values

#=========================================
def collectFinalStates(dr, dirs, dofName, timeStep):
  nDirs = len(dirs)
  data = np.zeros((1,nDirs))
  nDofs     = 0
  for i, iDir in enumerate(dirs):
    #fomFile = iDir + "/finalFomState_" + dofName + "_0"
    fomFile = iDir + "/state_timestep_" + str(timeStep) + "_" + dofName
    print("Reading {}".format(fomFile))
    if not dr:
      if not os.path.exists(fomFile):
        print("file {} does not exist".format(fomFile))
        sys.exit(1)

      fomData = np.loadtxt(fomFile, skiprows=1)
      nDofs = len(fomData)
      if len(data) == 1: data.resize((nDofs,nDirs))
      data[:,i] = fomData
  return data

#=========================================
def computeErrors(fomState, fomStateReconstructed, key, dofName):
  error = fomState-fomStateReconstructed
  print(key + "_{}_fom_minmax = {} {}".format(dofName, np.min(fomState), np.max(fomState)))
  print(key + "_{}_rom_minmax = {} {}".format(dofName, np.min(fomStateReconstructed), np.max(fomStateReconstructed)))
  print(key + "_{}_err_minmax = {} {}".format(dofName, np.min(error), np.max(error)))

  fomL2Norm, fomLinfNorm = la.norm(fomState), la.norm(fomState, np.inf)
  errL2Norm   = [ la.norm(error),         la.norm(error)/fomL2Norm ]
  errLinfNorm = [ la.norm(error, np.inf), la.norm(error, np.inf)/fomLinfNorm ]
  print(key + "_{}_err_abs_rel_ltwo_norms = {} {}".format(dofName, errL2Norm[0],   errL2Norm[1]))
  print(key + "_{}_err_abs_rel_linf_norms = {} {}".format(dofName, errLinfNorm[0], errLinfNorm[1]))


#=========================================
class Interpolator:
  def __init__(self, tag, nPts, nTests):
    self.supportedTags_ = ["nn", "linear", "quadratic", "cubicSpline", "splineK4", "rbf"]
    assert(tag in self.supportedTags_)
    self.tag_  = tag
    self.data_ = np.zeros((nPts, nTests))

  def getTag(self):
    return self.tag_

  def writeDataToFile(self, destDir, dofName):
    np.savetxt(destDir+'/'+self.tag_+'_'+dofName+'.txt', self.data_, fmt='%.18f')

  def viewData(self):
    return self.data_

  def evaluate(self, x, y, xtest, i):
    if (self.tag_ == "nn"):
      f=interp1d(x, y, kind='nearest', fill_value="extrapolate")

    elif (self.tag_ == "linear"):
      f=interp1d(x, y, kind='linear', fill_value="extrapolate")

    elif (self.tag_ == "quadratic"):
      f=interp1d(x, y, kind='quadratic', fill_value="extrapolate")

    elif (self.tag_ == "cubicSpline"):
      f=interp1d(x, y, kind='cubic', fill_value="extrapolate")

    elif (self.tag_ == "rbf"):
      f=Rbf(x, y)

    elif (self.tag_ == "splineK4"):
      # spline needs data in increasing order
      p = x.argsort(); x2 = x[p]; y2 = y[p]
      f = UnivariateSpline(x2, y2, k=4, ext=0, check_finite=True)

    self.data_[i][:] = f(xtest)


#=========================================
def doInterpolation(dr, scenario, workDir, dofName, trainDirs, testDirs, timeStep):
  # how many train dirs I am using
  numTrainDirs = len(trainDirs)

  # get train values used
  xTrain = getTargetValues(dr, args.scenario, trainDirs)
  print("\nTraining points {}".format(xTrain))
  assert(len(xTrain) == numTrainDirs)

  # get test values used
  xTest = getTargetValues(dr, args.scenario, testDirs)
  print("\nTest points {}".format(xTest))

  # collect train data: rows index the grid points, cols the training
  # so row (i,j) contains the value of dofName at the i-th grid point
  # obtained from the j-th training run
  trainData = collectFinalStates(dr, trainDirs, dofName, timeStep)
  nPts = trainData.shape[0]
  print(nPts)

  # create all interpolator objects
  myInterpolators = []
  myInterpolators.append(Interpolator("nn", nPts, len(xTest)))
  myInterpolators.append(Interpolator("linear", nPts, len(xTest)))
  if numTrainDirs >= 3:
    myInterpolators.append(Interpolator("quadratic", nPts, len(xTest)))
  if numTrainDirs >= 4:
    myInterpolators.append(Interpolator("cubicSpline", nPts, len(xTest)))
    myInterpolators.append(Interpolator("rbf", nPts, len(xTest)))
  if numTrainDirs >= 6: myInterpolators.append(Interpolator("splineK4", nPts, len(xTest)))

  # do actual interpolation for each mesh point
  print("\nInterpolating")
  for i in range(0, nPts):
    if i==int(nPts/4) or i==int(nPts/2) or i==int(3*nPts/4) or i==nPts-1:
      print("{} % done".format(float(i/nPts)*100))

    # for each row of the data, run interpolators
    for interp in myInterpolators: interp.evaluate(xTrain, trainData[i][:], xTest, i)

  # empty the train data since we dont need it anymore
  trainData = []

  # create folder for the interpolation
  interpDirFullPath = workDir+'/interpolation_n'+str(numTrainDirs)
  print("\nSaving interpolation results in {}".format(interpDirFullPath))
  if not os.path.exists(interpDirFullPath):
    os.system('mkdir -p ' + interpDirFullPath)

  # write interpolated data to file for each interpolator
  for interp in myInterpolators: interp.writeDataToFile(interpDirFullPath, dofName)

  # collect test data: rows index the grid points, cols the test runs
  # so row (i,j) contains the value of dofName at the i-th grid point
  # obtained from the j-th test run
  print("\nReading test data")
  testData = collectFinalStates(dr, testDirs, dofName, timeStep)
  assert(testData.shape[1] == len(xTest))

  # loop over tests and compute error
  for j in range(testData.shape[1]):
    print("\n")
    # compute errors
    for interp in myInterpolators:
      approxData = interp.viewData()[:,j]
      computeErrors(testData[:,j], approxData, interp.getTag()+"_test"+str(j), dofName)


###############################
if __name__== "__main__":
###############################
  parser = ArgumentParser()
  parser.add_argument("-dryrun", "--dryrun", "-dr", "--dr",
                      dest="dryRun", type=str2bool, default=True,
                      help="True: creates directory structures/files, does not run. Default=True.")

  parser.add_argument("-working-dir", "--working-dir", "-wdir", "--wdir",
                      dest="workDir", default="empty",
                      help="Target dir where to work. Must be set.")

  parser.add_argument("-dof-name", "--dof-name", "-dof", "--dof",
                      dest="dofName", default="empty",
                      help="Which dof: vp/sp to compute error for, must set.")

  parser.add_argument("-scenario", "--scenario",
                      dest="scenario", default=0, type=int,
                      help="Choices: 1 (uncertain shear velocity) or 2 (uncertain forcing period). Must set.")

  parser.add_argument("-ntrain", "--ntrain", "-numtrain", "--numtrain",
                      dest="nTrain", default=-1, type=int,
                      help="Num of training points to use. Must set.")

  parser.add_argument("-timestep", "--timestep",
                      dest="timeStep", default=-1, type=int,
                      help="Time step to use. Must set.")

  #------------------------------
  # parse all args
  #------------------------------
  args = parser.parse_args()
  assert(args.workDir != "empty")
  assert(args.dofName != "empty")
  assert(args.scenario == 1 or args.scenario == 2)
  assert(args.nTrain != -1)
  assert(args.timeStep != -1)

  # assert workDir exists
  assert(os.path.exists(args.workDir))

  # find all train dirs (sorted by run ID)
  trainDirs = getListOfDirs(args.dryRun, args.scenario, args.workDir, key="train", sortByRunID=True)
  # subselect only nTrain based on arg
  trainDirs = trainDirs[0:args.nTrain]
  print("\nTraining dirs:")
  for i in trainDirs: print(i)

  # find all FOM test dirs (sorted by run ID)
  testDirs = getListOfDirs(args.dryRun, args.scenario, args.workDir, key="test", sortByRunID=True)
  print("\nTest dirs:")
  for i in testDirs: print(i)

  # do interpolation
  doInterpolation(args.dryRun, args.scenario, args.workDir, args.dofName, trainDirs, testDirs, int(args.timeStep))
