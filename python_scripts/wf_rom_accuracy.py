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
from log_file_extractor import extractNumSteps

np.set_printoptions(edgeitems=10, linewidth=100000)

# here we only care about accuracy, so does not matter number of threads
# num of threads to use for fom
fomNThread = 36
romNThread = 18

#=========================================
def setSamplingFrequencyForTestRun(thisRunDic):
  # for testing, fix sampling such that we only collect 4 states
  finalRunTime=thisRunDic['general']['finalTime']
  timeStepSize=thisRunDic['general']['dt']
  totalSteps = finalRunTime/timeStepSize
  thisRunDic['io']['snapshotMatrix']['velocity']['freq'] = int(totalSteps/4)
  thisRunDic['io']['snapshotMatrix']['stress']['freq']   = int(totalSteps/4)

#=========================================
def createBaseDic(scenario, meshFullPath):
  dic = {
    'general':
    {'meshDir': meshFullPath,
     'dt': .1,
     'finalTime': 2000,
     'checkNumericalDispersion': True,
     'checkCfl': True,
     'includeMatPropInJacobian': True,
     'exploitForcingSparsity': True
    },
    # IO section
    'io':
    {
     'snapshotMatrix':
     {'binary': True,
      'velocity':  {'freq': 20, 'fileName': 'snaps_vp'},
      'stress':    {'freq': 20, 'fileName': 'snaps_sp'}},
      #
     'seismogram':
     {'binary': True, 'freq': 20, 'receivers': np.arange(5, 180, 5).tolist()}
    }
  }

  if scenario == 1:
    # recall that for scenario 1 we sample the velocity in second layer (see -1 below)
    dic['source'] = {
      'signal': {'kind': 'gaussDer', 'depth': 700., 'angle': 0., 'period': 45., 'delay': 140.}
    }
    dic['material'] = {
      'kind': 'bilayer',
      'layer1': {'density': [4500., 0.], 'velocity': [6000., 0.]},
      'layer2': {'depth': 700, 'density': [4500., 0.], 'velocity': [-1., 0.]}
    }

  elif scenario==2:
    # recall that for scenario 2 we sample the source signal period (see -1 below)
    dic['source'] = {
      'signal': {'kind': 'ricker', 'depth': 150., 'angle': 0., 'period': -1., 'delay': 90.}
    }
    dic['material'] = {'kind': 'prem'}

  return dic

#=========================================
def getSimulationTotalNumSteps(dic):
  finalRunTime=dic['general']['finalTime']
  timeStepSize=dic['general']['dt']
  return int(finalRunTime/timeStepSize)

#=========================================
def getRunIDFromDirName(dir):
  return int(dir.split('_')[-1])

#=========================================
def generateMesh(workDir, m):
  meshDir = workDir + "/meshes"
  if not os.path.exists(meshDir): os.system('mkdir ' + meshDir)
  generateMeshIfNeeded(meshDir, m[0], m[1])
  return meshDir+"/mesh"+str(m[0])+"x"+str(m[1])

#=========================================
def doFomRuns(dryRun, scenario, workDir, baseDic, exeDir, values, kind, runIDStart, reduceSampling=False):
  owd = os.getcwd()

  assert(kind=="train" or kind=="test")
  assert(runIDStart >= 0)
  print("values ", values)

  for i,it in enumerate(values):
    # clone yaml file from base case
    thisRunDic = copy.deepcopy(baseDic)

    if scenario==1:
      # replace velocity sample for this run
      thisRunDic['material']['layer2']['velocity'] = [float(it), 0.0]
    elif scenario==2:
      # replace the forcing period for this run
      thisRunDic['source']['signal']['period'] = float(it)

    # reduce the sampling freq for convenience if needed
    if reduceSampling==True: setSamplingFrequencyForTestRun(thisRunDic)

    # create dir for this run (append string to dirname to mark train/test and runID)
    runDirFullPath = createFomRunDirectory(workDir, fomNThread, thisRunDic, "_"+kind+"_"+str(i+runIDStart))
    if not os.path.exists(runDirFullPath): os.mkdir(runDirFullPath)

    # enter rundir, print yaml file and run
    os.chdir(runDirFullPath)
    with open('input.yaml', 'w') as yaml_file:
      yaml.dump(thisRunDic, yaml_file,  default_flow_style=False, sort_keys=False)

    if not dryRun:
      # run
      runExe(runDirFullPath, exeDir, fomExeName, fomNThread)

      # extract from the snapshots some target states that will be used for error checking
      linkExe(runDirFullPath, exeDir, extractStateExeName)
      samplingFreq = {'vp': thisRunDic['io']['snapshotMatrix']['velocity']['freq'],
                      'sp': thisRunDic['io']['snapshotMatrix']['stress']['freq']}

      totalTimeSteps = getSimulationTotalNumSteps(thisRunDic)

      for dof in ['vp', 'sp']:
        args = "./"+extractStateExeName +\
          " --snaps=./snaps_"+dof+"_0 binary" +\
          " --fsize="+str(1) +\
          " --outformat=ascii" +\
          " --timesteps=" + str(int(totalTimeSteps))+\
          " --samplingfreq=" + str(samplingFreq[dof]) +\
          " --outfileappend=" + dof

        print("Running: {} inside {}".format(args, runDirFullPath))
        logfile = 'out_extract_'+dof+'.log'
        os.system(args + ' > ' + logfile)
        # check that it worked
        assert('success' in open(logfile).read().split())
    else:
      print("dryrun:fom: {}".format(runDirFullPath))

    # return
    os.chdir(owd)

#=========================================
def runSVD(dryRun, numPodTrainPts, workDir, exeDir, binaryFlag):
  owd = os.getcwd()

  # find dirs with "train" in the name to use as trainings
  fomDirsFullPath = [workDir+'/'+d for d in os.listdir(workDir) if 'train' in d and 'fom' in d]
  # sort based on the ID
  def func(elem): return int(elem.split('_')[-1])
  fomDirsFullPath = sorted(fomDirsFullPath,key=func)
  fomDirsFullPath  = fomDirsFullPath[0:numPodTrainPts]

  # create the svd directory
  podDirFullPath = workDir+'/pod'
  if not os.path.exists(podDirFullPath): os.mkdir(podDirFullPath)

  # which svd method: recall that 0 is the old one, 1 usest A^T A
  svdMethod = '0'
  if (len(fomDirsFullPath) == 1):
    args = ("./"+svdExeName, workDir+'/'+fomDirsFullPath[0], binaryFlag, svdMethod)
  elif (len(fomDirsFullPath) == 2):
    args = ("./"+svdExeName, fomDirsFullPath[0], fomDirsFullPath[1], binaryFlag, svdMethod)
  elif (len(fomDirsFullPath) == 3):
    args = ("./"+svdExeName, fomDirsFullPath[0], fomDirsFullPath[1], fomDirsFullPath[2], binaryFlag, svdMethod)
  else:
    sys.exit("More than 3 train dirs not supported for svd")

  if not dryRun:
    # cd into dir where we do svd
    os.chdir(podDirFullPath)

    # link exeName from exeDir to workDir
    linkExe(podDirFullPath, exeDir, svdExeName)

    print("Running: {} inside {}".format(args, podDirFullPath))
    svdLogFilePath = podDirFullPath+'/out.log'

    logfile = open('out.log', 'w')
    p = subprocess.Popen(args, stdout=logfile, stderr=logfile)
    p.wait()
    logfile.close()

    # make sure that the svd did not fail by checking rank
    rankVp = extractRankFromSvdLogFile(svdLogFilePath, "vp")
    rankSp = extractRankFromSvdLogFile(svdLogFilePath, "sp")
    assert(rankVp!=0 and rankSp!=0)
    os.chdir(owd)
  else:
    print("dryrun:svd: {}".format(args))

  return podDirFullPath


#=========================================
def doRankOneRoms(dryRun, scenario, workDir, baseDic, exeDir, samples, podDirPath, modes):
  owd = os.getcwd()

  # where to find the basis
  vpLsv = podDirPath + '/lsv_vp'
  spLsv = podDirPath + '/lsv_sp'
  if not dryRun:
    assert( os.path.exists(vpLsv) )
    assert( os.path.exists(spLsv) )

  for i,it in enumerate(samples):
    # clone dic
    thisRunDic = copy.deepcopy(baseDic)

    if scenario==1:
      # then the values passed in refer to values of the velocity so replace in dir
      thisRunDic['material']['layer2']['velocity'] = [float(it), 0.0]
    elif scenario==2:
      # replace period ( delay is same as baseDic)
      thisRunDic['source']['signal']['period'] = float(it)

    # remove the seismogram section since for roms we do not compute it online
    seismoKey = 'seismogram'
    if seismoKey in thisRunDic['io']: del thisRunDic['io'][seismoKey]

    # set sampling freq for testing
    setSamplingFrequencyForTestRun(thisRunDic)

    # loop over modes
    for nModesVp, nModesSp in zip(modes['vp'], modes['sp']):
      # update the rom section
      thisRunDic['rom'] = {'velocity': {'numModes': int(nModesVp), 'modesFile': vpLsv},
                           'stress':   {'numModes': int(nModesSp), 'modesFile': spLsv}}

      # create dir for this run
      runDirFullPath = createRomRunDirectory(workDir, romNThread, thisRunDic, "_test_"+str(i))
      if not os.path.exists(runDirFullPath):
        os.mkdir(runDirFullPath)

      #enter, print yaml file and run
      os.chdir(runDirFullPath)
      with open('input.yaml', 'w') as yaml_file:
        yaml.dump(thisRunDic, yaml_file, default_flow_style=False, sort_keys=False)

      if not dryRun:
        # all these ROM runs are using rank-1 forcing, so force OpenMP affinity to spread
        runExe(runDirFullPath, exeDir, romExeName, romNThread, forceSpread=True)

        # reconstruct the FOM states that will be used for error checking
        linkExe(runDirFullPath, exeDir, reconstructFomStateExeName)
        samplingFreq = {'vp': thisRunDic['io']['snapshotMatrix']['velocity']['freq'],
                        'sp': thisRunDic['io']['snapshotMatrix']['stress']['freq']}
        totalTimeSteps = getSimulationTotalNumSteps(thisRunDic)

        for dof in ['vp', 'sp']:
          podModes = vpLsv    if dof=='vp' else spLsv
          nModes   = nModesVp if dof=='vp' else nModesSp
          args = "./"+reconstructFomStateExeName +\
            " --podmodes="+podModes+' binary' +\
            " --romsize="+str(nModes) +\
            " --romsnaps=./snaps_"+dof+"_0 binary" +\
            " --fsize="+str(1) +\
            " --outformat=ascii" +\
            " --timesteps=" + str(int(totalTimeSteps))+\
            " --samplingfreq=" + str(samplingFreq[dof]) +\
            " --outfileappend=" + dof

          print("Running: {} inside {}".format(args, runDirFullPath))
          logfile = 'out_reconstruct_'+dof+'.log'
          os.system(args + ' > ' + logfile)
          # check that it worked
          assert('success' in open(logfile).read().split())
      else:
        print("dryrun:rom: {}".format(runDirFullPath))

      # return
      os.chdir(owd)

#=========================================
def doRankTwoRoms(dryRun, scenario, workDir, baseDic, exeDir, samples, podDirPath, modes):
  owd = os.getcwd()

  # this is only ok for scenario 2 where we sample the forcing period
  assert(scenario==2)

  # where to find the basis
  vpLsv = podDirPath + '/lsv_vp'
  spLsv = podDirPath + '/lsv_sp'
  if not dryRun:
    assert( os.path.exists(vpLsv) )
    assert( os.path.exists(spLsv) )

  # loop over modes
  for nModesVp, nModesSp in zip(modes['vp'], modes['sp']):
    thisRunDic = copy.deepcopy(baseDic)

    # remove the seismogram section since for roms we do not compute it online
    if 'seismogram' in thisRunDic['io']: del thisRunDic['io']['seismogram']

    # the delay remains the same as set in baseDic, just update the periods
    # set sampling freq for testing
    setSamplingFrequencyForTestRun(thisRunDic)

    # replace periods so that we do a multiforcing run
    thisRunDic['sampling'] = {'params': ['signalPeriod'],
                              'values': samples.tolist(),
                              'forcingSize': int(len(samples))}

    # update the rom section
    thisRunDic['rom'] = {'velocity': {'numModes': int(nModesVp), 'modesFile': vpLsv},
                         'stress':   {'numModes': int(nModesSp), 'modesFile': spLsv}}

    # create dir for this run
    runDirFullPath = createRomRunDirectory(workDir, romNThread, thisRunDic)
    if not os.path.exists(runDirFullPath):os.mkdir(runDirFullPath)

    #enter, print yaml file and run
    os.chdir(runDirFullPath)
    with open('input.yaml', 'w') as yaml_file:
      yaml.dump(thisRunDic, yaml_file, default_flow_style=False, sort_keys=False)

    if not dryRun:
      runExe(runDirFullPath, exeDir, romExeName, romNThread)

      # reconstruct the FOM states that will be used for error checking
      linkExe(runDirFullPath, exeDir, reconstructFomStateExeName)
      samplingFreq = {'vp': thisRunDic['io']['snapshotMatrix']['velocity']['freq'],
                      'sp': thisRunDic['io']['snapshotMatrix']['stress']['freq']}
      totalTimeSteps = getSimulationTotalNumSteps(thisRunDic)
      fSz = thisRunDic['sampling']['forcingSize']

      for dof in ['vp', 'sp']:
        podModes = vpLsv    if dof=='vp' else spLsv
        nModes   = nModesVp if dof=='vp' else nModesSp
        args = "./"+reconstructFomStateExeName +\
          " --podmodes="+podModes+' binary' +\
          " --romsize="+str(nModes) +\
          " --romsnaps=./snaps_"+dof+"_0 binary" +\
          " --fsize="+str(fSz) +\
          " --outformat=ascii" +\
          " --timesteps=" + str(int(totalTimeSteps))+\
          " --samplingfreq=" + str(samplingFreq[dof]) +\
          " --outfileappend=" + dof

        print("Running: {} inside {}".format(args, runDirFullPath))
        logfile = 'out_reconstruct_'+dof+'.log'
        os.system(args + ' > ' + logfile)
        # check that it worked
        assert('success' in open(logfile).read().split())
    else:
      print("dryrun:rom: {}".format(runDirFullPath))

    # return
    os.chdir(owd)

#=========================================
def getListOfModesToRun(dryRun, podDirFullPath, energies):
  modes = {'vp': [], 'sp': []}
  for pct in energies:
    if (pct == 100.):
      svdLogFilePath = podDirFullPath+'/out.log'
      rankVp = extractRankFromSvdLogFile(svdLogFilePath, "vp") if not dryRun else 1
      rankSp = extractRankFromSvdLogFile(svdLogFilePath, "sp") if not dryRun else 1
      modes['vp'].append(np.max([rankVp, rankSp]))
      modes['sp'].append(np.max([rankVp, rankSp]))

    else:
      # run cumulative energy to compute various # of modes
      args   = ("python", "cumulative_energy.py",
                "-working-dir", podDirFullPath,
                "-percent", str(pct))
      popen  = subprocess.Popen(args, stdout=subprocess.PIPE); popen.wait()
      output = popen.stdout.read(); print("\n", output)
      # find vp
      res1 = re.search(cumEnRegExpVp, str(output))
      numBasisVp = int(res1.group().split()[1])
      assert(numBasisVp>0)
      # find sp
      res2 = re.search(cumEnRegExpSp, str(output))
      numBasisSp = int(res2.group().split()[1])
      assert(numBasisSp>0)

      #res = re.search(cumEnRegExp, str(output))
      #numBasis = int(res.group().split()[1])
      modes['vp'].append(numBasisVp)
      modes['sp'].append(numBasisSp)
  return modes

#=========================================
def computeRankOneRomErrors(dr, scenario, workDir):
  # in this case we have one dir for each test run so we should have
  # a matching number of fom and rom test runs where each run directory
  # has an test_id identifier

  # get all rom dirs
  romDirsFullPath = [workDir+'/'+d for d in os.listdir(workDir) if 'rom' in d]
  # get all fom test dirs
  fomDirsFullPath = [workDir+'/'+d for d in os.listdir(workDir) if 'fom' in d and 'test' in d]

  # loop over rom dirs
  for romDir in romDirsFullPath:
    # extract the test ID of this run
    id = getRunIDFromDirName(romDir)

    # find the corresponding fom test run
    reg = re.compile(r'.+fom_.+_test_'+str(id)+'$')
    mylist = list(filter(reg.match, fomDirsFullPath))
    assert(len(mylist)==1)
    fomDirFullPath = mylist[0]

    # extract total num of steps
    finalTimeStep = 1111 # dummy value
    if not dr:
      romLogFile = romDir+'/out.log'
      fomLogFile = fomDirFullPath+'/out.log'
      romFinalTimeStep = extractNumSteps(romLogFile)
      fomFinalTimeStep = extractNumSteps(fomLogFile)
      assert(romFinalTimeStep == fomFinalTimeStep)
      finalTimeStep = romFinalTimeStep

    # compute error at the final step
    for dof in ["vp", "sp"]:
      fomFile = fomDirFullPath+"/state_timestep_" + str(finalTimeStep) + "_" + dof
      romFile = romDir+"/fomReconstructedState_timestep_" + str(finalTimeStep) + "_" + dof

      owd = os.getcwd()
      args   = ("python", "compute_state_error.py",
                "-dryrun", str(dr),
                "-fomstate", fomFile,
                "-approxstate", romFile,
                "-dofname", dof)

      if not dr:
        print("Running: {}".format(args))
        logFileFullPath = romDir+'/finalError_'+dof+'.txt'
        logfile = open(logFileFullPath, 'w')
        p  = subprocess.Popen(args, stdout=logfile, stderr=logfile)
        p.wait()
        logfile.close()
      else:
        print("dryrun:rom_err: {}".format(args))
      os.chdir(owd)

#=========================================
def computeRankTwoRomErrors(dr, scenario, workDir):
  assert(scenario==2)

  # in this case we have fom tests enumrated as test_0, test_1, etc..
  # but we have one ROM dir for each # of modes that contains all samples of the forcing
  # inside each rom test dir, we have one file for each test as:
  # for instante, for numModes = X, we have a dir:
  #     rom_nThreads_..._nPod_X_X
  # inside we have:
  # finalFomState_{vp,sp}_i: final fom state for i-th test

  # get all fom test dirs
  fomDirsFullPath = [workDir+'/'+d for d in os.listdir(workDir) if 'fom' in d and 'test' in d]
  # get all rom dirs
  romDirsFullPath = [workDir+'/'+d for d in os.listdir(workDir) if 'rom' in d]

  # loop over fom test dirs
  for fomTestDir in fomDirsFullPath:
    # extract the test ID of this run
    testID = getRunIDFromDirName(fomTestDir)

    # loop over all rom dirs
    for romDir in romDirsFullPath:
      finalTimeStep = 1111 # dummy value
      if not dr:
        # extract total num of steps
        romLogFile = romDir+'/out.log'
        fomLogFile = fomTestDir+'/out.log'
        romFinalTimeStep = extractNumSteps(romLogFile)
        fomFinalTimeStep = extractNumSteps(fomLogFile)
        assert(romFinalTimeStep == fomFinalTimeStep)
        finalTimeStep = romFinalTimeStep

      # compute error
      for dof in ["vp", "sp"]:
        # inside the fom run I should have a single finalFomState_0 because this is a rank-one run
        fomFile = fomTestDir+"/state_timestep_" + str(finalTimeStep) + "_" + dof
        # inside the rom run I need to pick the right file because of the test ID
        romFile = romDir+"/fomReconstructedState_timestep_" + str(finalTimeStep) + "_f_" + str(testID) + "_" + dof

        owd = os.getcwd()
        args   = ("python", "compute_state_error.py",
                  "-dryrun", str(dr),
                  "-fomstate", fomFile,
                  "-approxstate", romFile,
                  "-dofname", dof)

        if not dr:
          print("Running: {}".format(args))
          logFileFullPath = romDir+'/finalError_'+dof+'_'+str(testID)+'.txt'
          logfile = open(logFileFullPath, 'w')
          p  = subprocess.Popen(args, stdout=logfile, stderr=logfile)
          p.wait()
          logfile.close()
        else:
          print("dryrun:rom_err: {}".format(args))
        os.chdir(owd)


#=========================================
def computeInterpolationErrors(dr, scenario, numPodTrainPts, workDir, finalTimeStep):
  # find out how many train dirs we have inside workDir
  trainDirs = [d for d in os.listdir(workDir) if "train" in d and "fom" in d]
  numTrains = len(trainDirs)
  assert(numTrains != 0)

  # check interpolation error for varying number of training points
  # even if it would be fair to only interpolate using the same data
  # used to generate the POD modes, we run interpolation for a few points
  # more to see how interpolation behaves
  for nTr in range(numPodTrainPts, numPodTrainPts+2, 1):
    print("Interpolation for nTr = {}".format(nTr))
    for dof in ["vp", "sp"]:
      owd = os.getcwd()
      args   = ("python", 'interpolate_state_target_step.py',
                "-dr", str(False),
                "-wdir", workDir,
                "-scenario", str(scenario),
                "-dof", dof,
                "-ntrain", str(nTr),
                "-timestep", str(finalTimeStep))

      if not dr:
        # the directory where the interpolation script puts data
        # fix this: if the interpolation script changes where interp data are stored,
        # this fails. So we need a more robust way to set the destination dir for the itnerpolation data
        interpDataDir = workDir+'/interpolation_n'+str(nTr)

        logFileFullPath = workDir+'/interpolation_'+dof+'.log'
        logfile = open(logFileFullPath, 'w')
        p  = subprocess.Popen(args, stdout=logfile, stderr=logfile)
        p.wait()
        logfile.close()

        if not os.path.exists(interpDataDir):
          sys.exit("The directory {} does not exist".format(interpDataDir))

        # move log file to interpDataDir
        os.system("mv " + logFileFullPath + ' ' + interpDataDir)
      else:
        print("dryrun:interp_err: {}".format(args))
        os.chdir(owd)

#========================================================
def createTrainTestSets(scenario, numPodTrainPts, pRange, totNumSamples, extrapPts):
  np.random.seed(6311276)

  # generate samples
  samples = np.random.uniform(low=float(pRange[0]), high=float(pRange[1]), size=(numSamples))
  print("Sampling values = ", samples)

  #----------------------
  # --- training set ---
  #----------------------
  # the values at the bounds are always the firt two points
  pTrainValues = np.array(pRange)
  if numPodTrainPts==3:
    # also add middle point
    middleValue = (pRange[0]+pRange[1])*0.5
    pTrainValues = np.append(pTrainValues,middleValue)

  # add random points in the interval
  pTrainValues = np.append(pTrainValues, np.random.choice(samples, int(0.60*numSamples), replace=False))

  #----------------------
  # --- test set ---
  #----------------------
  # take remaining from samples
  pTestValues = [i for i in samples if i not in pTrainValues]
  # append the extrapolation points
  pTestValues = np.append(pTestValues, extrapPts)

  # also use as test points the samples at the bounds (which are used for checking reproductive scenario)
  pTestValues = np.append(pTestValues, pRange[0])
  pTestValues = np.append(pTestValues, pRange[1])
  if numPodTrainPts==3:
    pTestValues = np.append(pTestValues, middleValue)

  print("Train values = ", pTrainValues)
  print("Test  values = ", pTestValues)
  return [pTrainValues, pTestValues]


###############################
if __name__== "__main__":
###############################
  parser = ArgumentParser()
  parser.add_argument("-dryrun", "--dryrun", "-dr", "--dr",
                      dest="dryRun", type=str2bool, default=True,
                      help="True: creates directory structures/files, does not run. Default=True.")

  parser.add_argument("-workingdir", "--workingdir", "-wdir", "--wdir",
                      dest="workDir", default="empty",
                      help="Target dir where to work. Required!")

  parser.add_argument("-exedir", "--exedir",
                      dest="exeDir", default="empty",
                      help="Exes dir is where I find exes. Required!")

  parser.add_argument("-scenario", "--scenario",
                      dest="scenario", type=int,
                      help="Choices: 1 (uncertain shear velocity), 2 (uncertain forcing period, fixed delay). Required!")

  parser.add_argument("-numsamples", "--numsamples",
                      dest="numSamples", type=int,
                      help="Num of samples. Required!")

  parser.add_argument("-numpodtrainpts", "--numpodtrainpts",
                      dest="numPodTrainPts", type=int,
                      help="Num of points for doing POD: current choices 2 (bounds), 3 (bounds and mid point). Required!")

  parser.add_argument("-mesh", "--mesh",
                      nargs="*", dest="meshSize", type=int,
                      help="Num of points to use along r and theta. Required!")

  #------------------------------
  # parse all args
  #------------------------------
  args = parser.parse_args()
  assert(args.workDir != "empty")
  assert(args.exeDir != "empty")
  assert(args.scenario in [1,2])
  assert(args.numSamples > 5)
  assert(args.numPodTrainPts in [2,3])
  workDir  = args.workDir
  scenario = args.scenario
  drun     = args.dryRun
  numPodTrainPts = args.numPodTrainPts
  assert( len(list(args.meshSize))==2 )

  # check if working dir exists, if not, make it
  if not os.path.exists(workDir): os.system('mkdir -p ' + workDir)

  # check that dir with exe exits
  if not os.path.exists(args.exeDir):
    sys.exit("The exedir {} provided is empty".format(args.exeDir))

  print("Running predictive for scenario {}".format(scenario))

  #------------------------------------------
  # list of energies for computing ROMs (100 % corresponds to the rank obtained from svd)
  energies = [99., 99.99, 99.999, 99.9999, 99.99999, 99.999999, 99.9999999, 99.99999999]

  # total samples to do (this is split over test/train)
  numSamples = args.numSamples

  # the parameter's range to use (depends on scenarion)
  pRange       = [6000., 6300.] if scenario==1 else [35., 65.]
  extrapValues = [5950., 6350.] if scenario==1 else [31., 69.]
  [pTrainValues, pTestValues] = createTrainTestSets(scenario, numPodTrainPts, pRange, numSamples, extrapValues)

  #------------------------------
  # running
  #------------------------------
  meshFullPath = generateMesh(workDir, args.meshSize)

  # the base dictionary to use
  baseDic = createBaseDic(scenario, meshFullPath)

  # -------------------------------------------
  # **** training runs for POD ****
  # -------------------------------------------
  # use runIDStart=0 so these runs will to be labeled as _train_0 and _train_1, ...
  print("\nDoing rom training runs")
  # make sure the bounds are always the first two points in the train set
  assert( pTrainValues[0] == pRange[0] and pTrainValues[1] == pRange[1] )
  doFomRuns(drun, scenario, workDir, baseDic, args.exeDir, pTrainValues[0:numPodTrainPts], kind="train", runIDStart=0)

  # -------------------------------------------
  # **** fom test runs ****
  # -------------------------------------------
  # use runIDStart=0 so these runs will to be labeled as _test_0, _test_1, _test_2, ...
  print("\nDoing fom test runs")
  doFomRuns(drun, scenario, workDir, baseDic, args.exeDir, pTestValues,
            kind="test", runIDStart=0, reduceSampling=True)

  # -------------------------------------------
  # **** run POD ****
  # -------------------------------------------
  print("\nDoing svd")
  binaryFlag     = '1' if baseDic['io']['snapshotMatrix']['binary'] else '0'
  #podDirFullPath = '/Users/fnrizzi/Desktop/waveSample/ex1/pod2'
  #podDirFullPath = '/home/fnrizzi/waveRuns/romAccuracyScenario2/train2/pod'
  podDirFullPath = runSVD(drun, numPodTrainPts, workDir, args.exeDir, binaryFlag)

  # compute list of modes corresponding to target energies
  # for dryRun, just put dummy values 1,2 to test the workflow but nothing is run anyway
  modes = getListOfModesToRun(drun, podDirFullPath, energies) if not drun else [1, 2]
  print('List of modes {}'.format(modes))

  # -------------------------------------------
  # **** test runs ROM ****
  # -------------------------------------------
  print("\nDoing rom test runs")
  if (scenario==1):
    #scenario 1 deals with shear velocity in layer2, so we need to do multiple roms, one per sample
    doRankOneRoms(drun, scenario, workDir, baseDic, args.exeDir, pTestValues, podDirFullPath, modes)
    print("\nComputing rom errors")
    computeRankOneRomErrors(drun, scenario, workDir)
  else:
    # for scenario 2 we can do a single rom run with multiforcing, gives same result as doing multiple runs
    doRankTwoRoms(drun, scenario, workDir, baseDic, args.exeDir, pTestValues, podDirFullPath, modes)
    print("\nComputing rom errors")
    computeRankTwoRomErrors(drun, scenario, workDir)

  # -------------------------------------------
  # **** other train runs for interpolation ****
  # -------------------------------------------
  print("\nDoing remaining training runs")
  doFomRuns(drun, scenario, workDir, baseDic, args.exeDir, pTrainValues[numPodTrainPts:], kind="train",
            runIDStart=numPodTrainPts, reduceSampling=True)

  # compute interpolation errors
  print("\nComputing interpolation errors")
  totalTimeSteps = getSimulationTotalNumSteps(baseDic)
  computeInterpolationErrors(drun, scenario, numPodTrainPts, workDir, totalTimeSteps)
