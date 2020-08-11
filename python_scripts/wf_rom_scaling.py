#!/usr/bin/env python

from random import randrange, uniform
import re, sys, os, time, yaml, copy
import numpy as np
from argparse import ArgumentParser
import shutil, subprocess
import matplotlib.pyplot as plt
from constants import *
from regexes import *
from utils import *
from create_run_directory import createFomRunDirectory, createRomRunDirectory
from log_file_extractor import *

#=========================================
def createBaseDic(meshFullPath):
  dic = {
    'general':
    {'meshDir': meshFullPath,
     'dt': .1, 'finalTime': 400.,
     'checkNumericalDispersion': True,
     'checkCfl': True,
     'includeMatPropInJacobian': True,
     'exploitForcingSparsity': True
    },
    # source and material (for rom scaling this is not really used since we use dummy values)
    'source':
    {'signal':
     {'kind': 'gaussDer', 'depth': 640., 'angle': 0., 'period': 55., 'delay': 170.}
    }
  }

  dic['material'] = {'kind': 'prem'}
  return dic

#=========================================
def generateMeshes(meshParentDir, m):
  if not os.path.exists(meshParentDir): os.system('mkdir ' + meshParentDir)
  generateMeshIfNeeded(meshParentDir, m[0], m[1])

#=========================================
def doRankOneRomScaling(dryRun, workDir, baseDic, exeDir, modes, threads):
  owd = os.getcwd()

  # loop over modes
  for i, nModes in enumerate(modes):
    # clone dictionary
    thisRunDic = copy.deepcopy(baseDic)

    #(dont provide svd files so it runs dummy rom)
    thisRunDic['rom'] = {'velocity': {'numModes': int(nModes)},
                         'stress':   {'numModes': int(nModes)}}

    # also turns off ROM jacobians
    thisRunDic['rom']['disableCompRomJacs'] = True;

    # loop over threads
    for nThread in threads:
      # create dir for this run
      runDirFullPath = createRomRunDirectory(workDir, nThread, thisRunDic)
      if not os.path.exists(runDirFullPath): os.mkdir(runDirFullPath)

      #enter, print yaml file and run
      os.chdir(runDirFullPath)
      with open('input.yaml', 'w') as yaml_file:
        yaml.dump(thisRunDic, yaml_file, default_flow_style=False, sort_keys=False)

      if not dryRun:
        # for rank-1 forcing, always force OpenMP affinity to spread
        runExe(runDirFullPath, exeDir, romExeName, nThread, forceSpread=True)
      else:
        print("dryrun:rom: {}".format(runDirFullPath))

      # return
      os.chdir(owd)


#=========================================
def doRankTwoRomScaling(dryRun, workDir, baseDic, exeDir, modes, threads, fSizes):
  owd = os.getcwd()

  # range for sampling the forcing period (this does not really matter since for scaling
  # we do not care about real values, but we need to set it anyway)
  pRange = [60., 70.]

  # fix seed for samples of forcing period
  np.random.seed(4345454)

  # loop over f sizes
  for fSize in fSizes:
    # list of periods: pick as many samples as  size so that we can handle all in one run
    fPeriods = np.random.uniform( low=float(pRange[0]), high=float(pRange[1]), size=fSize )

    # clone dic
    thisRunDic = copy.deepcopy(baseDic)

    # I need the add the sampling section to sample period
    thisRunDic['sampling'] = {'params': ['signalPeriod'],
                              'values': fPeriods.tolist(),
                              'forcingSize': int(fSize)}

    # set delay based on max period
    maxT = np.max(fPeriods)
    thisRunDic['source']['signal']['delay'] = float(maxT)*3.

    # loop over modes
    for i, nModes in enumerate(modes):
      #(dont prodive svd files so it runs dummy rom)
      thisRunDic['rom'] = {'velocity': {'numModes': int(nModes)},
                           'stress':   {'numModes': int(nModes)}}

      # also turns off ROM jacobians
      thisRunDic['rom']['disableCompRomJacs'] = True;

      # loop over threads
      for nThread in threads:
        # create dir for this run
        runDirFullPath = createRomRunDirectory(workDir, nThread, thisRunDic)
        if not os.path.exists(runDirFullPath): os.mkdir(runDirFullPath)

        #enter, print yaml file and run
        os.chdir(runDirFullPath)
        with open('input.yaml', 'w') as yaml_file:
          yaml.dump(thisRunDic, yaml_file, default_flow_style=False, sort_keys=False)

        if not dryRun:
          if (fSize==1):
            runExe(runDirFullPath, exeDir, romExeName, nThread, forceSpread=True)
          else:
            runExe(runDirFullPath, exeDir, romExeName, nThread)

          # assert that rankTwo was enabled
          assert(extractRankTwoFlag(runDirFullPath+'/out.log'))
        else:
          print("dryrun:rom: {}".format(runDirFullPath))

        # return
        os.chdir(owd)

#=====================================================================
def findRomRunDirectory(workDir, nthreads, fSize, nModes):
  reg = re.compile(r'.*_nThreads_'+str(nthreads)+'_.*_fRank_'+str(fSize)+'_nPod_'+str(nModes))
  dirs = [d for d in os.listdir(workDir)]
  newlist = list(filter(reg.match, dirs))
  print("New list", newlist)
  assert(len(newlist)==1)
  print("\n")
  return newlist[0]

#=====================================================================
def extractData(workDir, threads, fSizes, modes, outFileFullPath):
  # col0        = nthreads
  # col1        = fSize
  # col2        = num modes
  # col3        = # dofs
  # col4        = CI (flops/bytes)
  # col5,6,7    = [ave, min, max BW] (GB/s)
  # col8,9,10   = [ave, min, max ] (Gflops)
  # col11,12,13 = [ave, min, max time] (ms)
  # col14,15    = [loopTime, collectionTime] (sec)
  # col16       = [jacobTime] (sec)
  # col17       = [finalTime] (sec)
  # col18       = [numSteps]

  numThreads = len(threads)
  numFs      = len(fSizes)
  numModes   = len(modes)
  numRows = numFs * numThreads * numModes

  numCols = 19
  data = np.zeros((numRows, numCols))

  iR = 0
  for nM in modes:
    for nThread in threads:
      for fS in fSizes:
        print("extracting data for fSize {} with {} threads and {} modes".format(fS, nThread, nM))

        # find inside workDir the dir matching mesh and threads
        dirPath = findRomRunDirectory(workDir, nThread, fS, nM)
        print("found dir {} ".format(dirPath))

        # the log file
        fomLogFile = workDir+'/'+dirPath+'/out.log'
        print(fomLogFile)

        # nThreads
        data[iR][0] = nThread
        # fSize
        data[iR][1] = fS
        # nModes
        data[iR][2] = nM
        # get num dofs (vp + sp)
        data[iR][3] = np.int64(extractNumTotalDofs(fomLogFile))
        # get flops/bytes
        data[iR][4] = extractComputIntensity(fomLogFile)
        # get mem BW
        data[iR][5:8] = extractMemBw(fomLogFile)
        # get Gflops
        data[iR][8:11] = extractGflops(fomLogFile)
        # get perf times
        data[iR][11:14] = extractPerfTimes(fomLogFile)
        # get loop times
        data[iR][14] = extractLoopTime(fomLogFile)
        data[iR][15] = extractDataIoTime(fomLogFile)
        # get jacobian time
        data[iR][16] = extractRomJacobianTime(fomLogFile)
        # get final time
        data[iR][17] = extractFinalTime(fomLogFile)
        # get numSteps
        data[iR][18] = extractNumSteps(fomLogFile)

        iR += 1

  np.savetxt(outFileFullPath, data,
             fmt=['%d', '%d', '%d', '%d',
                  '%.8f',
                  '%9.4f','%9.4f','%9.4f',
                  '%9.4f','%9.4f','%9.4f',
                  '%6.3f','%6.3f','%6.3f',
                  '%9.4f','%9.4f','%9.4f',
                  '%9.4f','%9.4f'])


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
  if not os.path.exists(args.workDir): os.system('mkdir -p ' + args.workDir)

  # check that dir with exe exits
  if not os.path.exists(args.exeDir):
    sys.exit("The exedir {} provided is empty".format(args.exeDir))

  #------------------------------------------
  # this script tests the ROM scaling wrt threads, for given f and # modes:
  # - mesh is fixed (mesh has impact on the memory footprint and comp of romJacobians but not on rom itself)
  # - we turn off computation of the jacobians just set them to random to make things faster
  # - we use dummy basis (just random values) for pure scaling of the rom we do not care about accuracy
  # - IO is off here for convenience
  # - we run for a limited number of step since we do not care of having real data

  mesh         = [300, 1200]
  meshDir      = args.workDir + "/meshes"
  meshFullPath = meshDir+"/mesh"+str(mesh[0])+"x"+str(mesh[1])

  # the base dictionary to use
  baseDic = createBaseDic(meshFullPath)

  # list of threads case
  threads = [1, 2, 4, 8, 12, 18, 36, 72]

  # for doing rom scaling with rank-2 forcing, we need to set the sizes of the forcing
  fSizes  = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]

  # list of modes
  modes = [64, 128, 256, 512, 1024, 2048, 4096]
  print('List of modes {}'.format(modes))

  #------------------------------
  # running
  #------------------------------
  # generateMeshes(meshDir, mesh)

  # # rank-1 rom scaling
  # print("\nDoing rom rank-1 scaling runs")
  # doRankOneRomScaling(args.dryRun, args.workDir, baseDic, args.exeDir, modes, threads)

  # # rank-2 rom scaling
  # print("\nDoing rom rank-2 scaling runs")
  # doRankTwoRomScaling(args.dryRun, args.workDir, baseDic, args.exeDir,
  #                     modes, threads, fSizes)

  if not args.dryRun:
    # when parsing data, also parse the single forcing case
    fSizesToParse = [1] + fSizes
    extractData(args.workDir, threads, fSizesToParse, modes, args.workDir+'/scaling.txt')
