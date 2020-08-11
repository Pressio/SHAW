#!/usr/bin/env python

from random import randrange, uniform
import re, sys, os, time, yaml, copy
import numpy as np
from argparse import ArgumentParser
import shutil, subprocess
from constants import *
from utils import *
from create_run_directory import createFomRunDirectory, createRomRunDirectory
from log_file_extractor import *

#=========================================
def createBaseDic():
  baseDic = {
    'general' :
    {'meshDir': "empty", 'dt': 0.05, 'finalTime': 50.,
     'checkNumericalDispersion': True,
     'checkCfl': True,
     'includeMatPropInJacobian': True,
     'exploitForcingSparsity': True
    },
    'source':
    {'signal':
     {'kind': 'gaussDer', 'depth': 640., 'angle': 0., 'period': 60., 'delay': 180.}
    },
    'material':
    {'kind': 'unilayer',
     'layer': {'density': [2500., 0.], 'velocity': [5000., 0.]}
    }
  }
  return baseDic


#=====================================================================
def generateMeshes(meshParentDir, cases):
  if not os.path.exists(meshParentDir):
    os.system('mkdir ' + meshParentDir)
  for m in cases:
    generateMeshIfNeeded(meshParentDir, m[0], m[1])

#=====================================================================
def doRunsRankOneForcing(dryRun, workDir, baseDic, exeDir, meshParentDir, threads):
  # list all meshes
  meshes = [d for d in os.listdir(meshParentDir)]
  print(meshes)
  owd = os.getcwd()

  # loop over threads
  for nThread in threads:
    # loop over meshes
    for iMesh in meshes:
      fullMeshPath = meshParentDir+"/"+iMesh

      # clone yaml file from base case, replace mesh dir
      thisRunDic = copy.deepcopy(baseDic)
      thisRunDic['general']['meshDir'] = fullMeshPath

      # create dir for this run
      runDirFullPath = createFomRunDirectory(workDir, nThread, thisRunDic)
      if not os.path.exists(runDirFullPath):
        os.mkdir(runDirFullPath)

      #enter, print yaml file and run
      os.chdir(runDirFullPath)
      with open('input.yaml', 'w') as yaml_file:
        yaml.dump(thisRunDic, yaml_file,
                  default_flow_style=False, sort_keys=False)

      if not dryRun:
        runExe(runDirFullPath, exeDir, fomExeName, nThread)
      else:
        print("dryrun:fom: {}".format(runDirFullPath))

      # return
      os.chdir(owd)

#=====================================================================
def doRunsRankTwoForcing(dryRun, workDir, baseDic, exeDir, meshParentDir, threads, fSizes):
  # list all meshes
  meshes = [d for d in os.listdir(meshParentDir)]
  print(meshes)
  owd = os.getcwd()

  # loop over f sizes
  for fSize in fSizes:
    # list of periods: pick as many samples as f size so that we can handle all in one run
    fPeriods = np.linspace(60, 80, fSize)
    # mean value of period
    aveT = np.mean(fPeriods)

    # loop over threads
    for nThread in threads:
      # loop over meshes
      for iMesh in meshes:
        fullMeshPath = meshParentDir+"/"+iMesh

        # clone yaml file from base case
        thisRunDic = copy.deepcopy(baseDic)

        # replace mesh dir
        thisRunDic['general']['meshDir'] = fullMeshPath

        # I need the add the sampling section
        thisRunDic['sampling'] = {'params': ['signalPeriod'],
                                  'values': fPeriods.tolist(),
                                  'forcingSize': int(fSize)}
        # set delay based on mean period
        thisRunDic['source']['signal']['delay'] = float(aveT)*3.

        # create dir for this run
        runDirFullPath = createFomRunDirectory(workDir, nThread, thisRunDic)
        if not os.path.exists(runDirFullPath):
          os.mkdir(runDirFullPath)

        #enter, print yaml file and run
        os.chdir(runDirFullPath)
        with open('input.yaml', 'w') as yaml_file:
          yaml.dump(thisRunDic, yaml_file,
                    default_flow_style=False, sort_keys=False)

        if not dryRun:
          runExe(runDirFullPath, exeDir, fomExeName, nThread)
        else:
          print("dryrun:fom: {}".format(runDirFullPath))

        # return
        os.chdir(owd)


#=====================================================================
def findFomRunDirectory(workDir, nvp, nth, nthreads, fSize):
  reg = re.compile(r'.*mesh'+str(nvp)+'x'+str(nth)+'_nThreads_'+str(nthreads)+'.*_fRank_'+str(fSize)+'$')
  dirs = [d for d in os.listdir(workDir)]
  print("Full list", dirs)
  print("\n")
  newlist = list(filter(reg.match, dirs))
  print("New list", newlist)
  print("\n")
  return newlist[0]

#=====================================================================
def extractData(workDir, meshList, threads, dataFileFullPath, fSizes):
  # col0        = nthreads
  # col1        = fSize
  # col2        = # dofs
  # col3        = CI (flops/bytes)
  # col4,5,6    = [ave, min, max BW] (GB/s)
  # col7,8,9    = [ave, min, max ] (Gflops)
  # col10,11,12 = [ave, min, max time] (ms)
  # col13,14    = [loopTime, collectionTime] (sec)
  # col15       = [finalTime] (sec)
  # col16       = [numSteps]

  numMeshes  = len(meshList)
  numThreads = len(threads)
  numRows = numMeshes * numThreads * len(fSizes)
  numCols = 17
  data = np.zeros((numRows, numCols))

  iR = 0
  for iMesh in meshList[:]:
    for nThread in threads[:]:
      for fS in fSizes:
        print("extracting data for mesh {} with {} threads and {} fSize".format(iMesh, nThread, fS))

        # find inside workDir the dir matching mesh and threads
        fomDirPath = findFomRunDirectory(workDir, iMesh[0], iMesh[1], nThread, fS)
        print("found dir {} ".format(fomDirPath))

        # the log file
        fomLogFile = workDir+'/'+fomDirPath+'/out.log'
        print(fomLogFile)

        # nThreads
        data[iR][0] = nThread
        # fSize
        data[iR][1] = fS
        # get num dofs (vp + sp)
        data[iR][2] = np.int64(extractNumTotalDofs(fomLogFile))
        # get flops/bytes
        data[iR][3] = extractComputIntensity(fomLogFile)
        # get mem BW
        data[iR][4:7] = extractMemBw(fomLogFile)
        # get Gflops
        data[iR][7:10] = extractGflops(fomLogFile)
        # get perf times
        data[iR][10:13] = extractPerfTimes(fomLogFile)
        # get loop times
        data[iR][13] = extractLoopTime(fomLogFile)
        data[iR][14] = extractDataIoTime(fomLogFile)
        # get final time
        data[iR][15] = extractFinalTime(fomLogFile)
        # get numSteps
        data[iR][16] = extractNumSteps(fomLogFile)

        iR += 1

  np.savetxt(dataFileFullPath, data,
             fmt=['%d', '%d', '%d',
                  '%.8f',
                  '%9.4f','%9.4f','%9.4f',
                  '%9.4f','%9.4f','%9.4f',
                  '%6.3f','%6.3f','%6.3f',
                  '%9.4f', '%9.4f',
                  '%9.4f', '%9.4f'])

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

  parser.add_argument("-exe-dir", "--exe-dir", "-exedir", "--exedir",
                      dest="exeDir", default="empty",
                      help="Exes dir is where I find exes. Must be set.")

  #------------------------------
  # parse all args
  #------------------------------
  args = parser.parse_args()
  assert(args.workDir != "empty")
  assert(args.exeDir != "empty")

  # check if working folder exists, if not, make it
  if not os.path.exists(args.workDir):
    os.system('mkdir -p ' + args.workDir)

  if not os.path.exists(args.exeDir):
    sys.exit("The exedir {} provided is empty".format(args.exeDir))

  #------------------------------
  # [nr, ntheta]
  meshList = [ [256, 1024], [512, 2048], [1024, 4096], [2048, 8192] ]

  # list of threads case
  threads  = [2, 4, 8, 12, 18, 36, 72]

  # list of forcing sizes
  fSizes  = [2, 4, 8, 16, 32, 48]

  meshesDir = args.workDir + "/meshes"
  #generateMeshes(meshesDir, meshList)

  # the base dictionary to use
  baseDic = createBaseDic()

  #------------------------------
  #doRunsRankOneForcing(args.dryRun, args.workDir, baseDic, args.exeDir, meshesDir, threads)
  #doRunsRankTwoForcing(args.dryRun, args.workDir, baseDic, args.exeDir, meshesDir, threads, fSizes)

  #if not dr:
  extractData(args.workDir, meshList, threads, args.workDir+'/scaling.txt', [1]+fSizes)
