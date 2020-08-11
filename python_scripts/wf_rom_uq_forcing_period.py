#!/usr/bin/env python

from random import randrange, uniform
import re, sys, os, time, yaml, copy, glob
import numpy as np
from argparse import ArgumentParser
import shutil, subprocess
import matplotlib.pyplot as plt
from constants import *
from regexes import *
from utils import *
from create_run_directory import createFomRunDirectory, createRomRunDirectory
from log_file_extractor import extractNumSteps

from wf_rom_accuracy import generateMesh, runSVD, getListOfModesToRun

np.set_printoptions(edgeitems=10, linewidth=100000)

# here we only care about accuracy, so does not matter number of threads
# num of threads to use for fom
fomNThread = 4
romNThread = 4

#=========================================
def setSamplingFrequencyForTestRunFom(thisRunDic):
  # for testing, we don't care about the snapshot matrix,
  # we only care for seismogram so fix sampling such
  # that we only collect 2 states
  finalRunTime=thisRunDic['general']['finalTime']
  timeStepSize=thisRunDic['general']['dt']
  totalSteps = finalRunTime/timeStepSize
  thisRunDic['io']['snapshotMatrix']['velocity']['freq'] = int(totalSteps/2)
  thisRunDic['io']['snapshotMatrix']['stress']['freq']   = int(totalSteps/2)

#=========================================
def setSamplingFrequencyForTrainRunFom(thisRunDic):
  finalRunTime=thisRunDic['general']['finalTime']
  timeStepSize=thisRunDic['general']['dt']
  totalSteps = finalRunTime/timeStepSize
  # for train fom run we need to have a higher sampling freq
  # because we need to compute the POD modes
  thisRunDic['io']['snapshotMatrix']['velocity']['freq'] = int(20)
  thisRunDic['io']['snapshotMatrix']['stress']['freq']   = int(20)

#=========================================
def setSamplingFrequencyForRunRom(thisRunDic):
  finalRunTime=thisRunDic['general']['finalTime']
  timeStepSize=thisRunDic['general']['dt']
  totalSteps = finalRunTime/timeStepSize
  # in theory we should set it here equal to those used for the
  # seismogram of the fom so that we have an exact correspondence
  thisRunDic['io']['snapshotMatrix']['velocity']['freq'] = int(10)
  thisRunDic['io']['snapshotMatrix']['stress']['freq']   = int(10)

#=========================================
def createBaseDicProductionRun(meshFullPath):
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
     {'binary': False, 'freq': 10, 'receivers': [10., 30., 60., 90., 120., 150.]} #np.arange(5, 180, 5).tolist()}
    }
  }

  # recall that we sample the source signal period (see -1 below)
  dic['source'] = {
    'signal': {'kind': 'ricker', 'depth': 150., 'angle': 0., 'period': -1., 'delay': 90.}
  }
  dic['material'] = {'kind': 'prem'}

  return dic

#=========================================
def createBaseDicTestingOnly(meshFullPath):
  dic = {
    'general':
    {'meshDir': meshFullPath,
     'dt': .25,
     'finalTime': 2000,
     'checkNumericalDispersion': False,
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
     {'binary': False, 'freq': 10, 'receivers': [10., 30., 60., 90., 120., 150.]} #np.arange(5, 180, 5).tolist()}
    },
    # material
    'material': {'kind': 'prem'}
  }
  # recall that we sample the source signal period (see -1 below)
  dic['source'] = {
    'signal': {'kind': 'gaussDer', 'depth': 150., 'angle': 0., 'period': -1., 'delay': 210.}
  }
  return dic


#=========================================
def doFomRunsRank1(dryRun, workDir, baseDic, exeDir, values, trainOrTest, runIDStart):
  owd = os.getcwd()
  assert(trainOrTest=="train" or trainOrTest=="test")
  assert(runIDStart >= 0)
  print("values ", values)

  for i,it in enumerate(values):
    # clone yaml file from base case
    thisRunDic = copy.deepcopy(baseDic)

    # replace the forcing period for this run
    thisRunDic['source']['signal']['period'] = float(it)

    if trainOrTest=="train":  setSamplingFrequencyForTrainRunFom(thisRunDic)
    elif trainOrTest=="test": setSamplingFrequencyForTestRunFom(thisRunDic)

    # create dir for this run (append string to dirname to mark train/test and runID)
    runDirFullPath = createFomRunDirectory(workDir, fomNThread, thisRunDic,
                                           "_"+trainOrTest+"_"+str(i+runIDStart))
    if not os.path.exists(runDirFullPath): os.mkdir(runDirFullPath)

    # enter rundir, print yaml file and run
    os.chdir(runDirFullPath)
    with open('input.yaml', 'w') as yaml_file:
      yaml.dump(thisRunDic, yaml_file,  default_flow_style=False, sort_keys=False)

    if not dryRun:
      # run
      runExe(runDirFullPath, exeDir, fomExeName, fomNThread)
    else:
      print("dryrun:fom: {}".format(runDirFullPath))

    os.chdir(owd)

#=========================================
def doFomRunsRank2(dryRun, workDir, baseDic, exeDir, values,
                   forcingSizeForRank2Fom, trainOrTest, runIDStart):
  owd = os.getcwd()
  assert(trainOrTest=="train" or trainOrTest=="test")
  assert(runIDStart >= 0)
  print("values ", values)

  # clone yaml file from base case
  thisRunDic = copy.deepcopy(baseDic)

  if trainOrTest=="train":  setSamplingFrequencyForTrainRunFom(thisRunDic)
  elif trainOrTest=="test": setSamplingFrequencyForTestRunFom(thisRunDic)

  # the delay remains the same as set in baseDic, just update the periods
  # replace periods so that we do a multiforcing run
  assert( len(values) % forcingSizeForRank2Fom == 0)
  if len(values) >= forcingSizeForRank2Fom:
    actualForcingSize = forcingSizeForRank2Fom
  else:
    actualForcingSize = len(values)
  thisRunDic['sampling'] = {'params': ['signalPeriod'],
                            'values': values.tolist(),
                            'forcingSize': int(actualForcingSize)}

  # create dir for this run (append string to dirname to mark train/test and runID)
  runDirFullPath = createFomRunDirectory(workDir, fomNThread, thisRunDic,
                                         "_"+trainOrTest+"_"+str(runIDStart))
  if not os.path.exists(runDirFullPath): os.mkdir(runDirFullPath)

  # enter rundir, print yaml file and run
  os.chdir(runDirFullPath)
  with open('input.yaml', 'w') as yaml_file:
    yaml.dump(thisRunDic, yaml_file,  default_flow_style=False, sort_keys=False)

  if not dryRun:
    runExe(runDirFullPath, exeDir, fomExeName, fomNThread)
    # don't need stress data
    os.system('rm -rf ' + runDirFullPath+'/snaps_sp_*')
  else:
    print("dryrun:fom: {}".format(runDirFullPath))

  os.chdir(owd)


#=========================================
def doRankTwoRoms(dryRun, workDir, baseDic, exeDir, samples,
                  maxForcingSizeForRank2Rom, podDirPath, modes):
  owd = os.getcwd()
  # where to find the basis
  vpLsv = podDirPath + '/lsv_vp'
  spLsv = podDirPath + '/lsv_sp'
  if not dryRun:
    assert( os.path.exists(vpLsv) )
    assert( os.path.exists(spLsv) )

  # loop over modes
  for nModesVp, nModesSp in zip(modes['vp'], modes['sp']):
    thisRunDic = copy.deepcopy(baseDic)

    setSamplingFrequencyForRunRom(thisRunDic)

    # remove the seismogram section since for roms we do not compute it online
    if 'seismogram' in thisRunDic['io']: del thisRunDic['io']['seismogram']

    # the delay remains the same as set in baseDic, just update the periods
    # replace periods so that we do a multiforcing run
    #assert( len(samples) % maxForcingSizeForRank2Rom == 0)
    if len(samples) >= maxForcingSizeForRank2Rom:
      actualForcingSize = maxForcingSizeForRank2Rom
    else:
      actualForcingSize = len(samples)
    thisRunDic['sampling'] = {'params': ['signalPeriod'],
                              'values': samples.tolist(),
                              'forcingSize': int(actualForcingSize)}

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
      # don't need stress data
      os.system('rm -rf ' + runDirFullPath+'/snaps_sp_*')
    else:
      print("dryrun:rom: {}".format(runDirFullPath))

    # return
    os.chdir(owd)

#=========================================
def extractRomSizeFromInputFile(dr, dirPath):
  with open(dirPath+'/input.yaml') as file:
    # The FullLoader parameter handles the conversion from YAML
    # scalar values to Python the dictionary format
    ifile = yaml.load(file, Loader=yaml.FullLoader)
    return int( ifile['rom']['velocity']['numModes'] )

#=========================================
def extractForcingSizeFromInputFile(dr, dirPath):
  with open(dirPath+'/input.yaml') as file:
    # The FullLoader parameter handles the conversion from YAML
    # scalar values to Python the dictionary format
    ifile = yaml.load(file, Loader=yaml.FullLoader)
    return int( ifile['sampling']['forcingSize'] )

#=========================================
def extractSeismogramGridPoints(fomDir):
  # load the log file
  logFile = fomDir+'/out.log'

  reg = re.compile(r'Point mapped to GID = \d+')
  file1 = open(logFile, 'r')
  strings = re.findall(reg, file1.read())
  file1.close()
  assert(strings)
  print(strings)
  gids = [i.split()[5] for i in strings]
  return gids


#=========================================
def reconstructSeismogramRom(dryRun, workDir, exeDir, samples, podDirPath):
  owd = os.getcwd()
  # where to find the basis
  vpLsv = podDirPath + '/lsv_vp'
  spLsv = podDirPath + '/lsv_sp'
  if not dryRun:
    assert( os.path.exists(vpLsv) )
    assert( os.path.exists(spLsv) )

  # find FOM test dir, should be only one
  fomDirsFullPath = [workDir+'/'+d for d in os.listdir(workDir) if 'fom' in d and 'test' in d]
  assert(len(fomDirsFullPath)==1)
  print(fomDirsFullPath)

  # from the FOM test run, extract all grid points where the seismogram is computed
  gids = extractSeismogramGridPoints(fomDirsFullPath[0])
  gidsString = ' '.join(gids)
  print(gidsString)

  # find all rom dirs
  romDirsFullPath = [workDir+'/'+d for d in os.listdir(workDir) if 'rom' in d]
  # sort based on the rom size which should be the last item in the dir name
  def func(elem): return int(elem.split('_')[-1])
  romDirsFullPath = sorted(romDirsFullPath,key=func)
  print(romDirsFullPath)

  # loop over each ROM dir
  for romDir in romDirsFullPath:
    #enter dir
    os.chdir(romDir)

    # extract the rom size
    romSize = extractRomSizeFromInputFile(dryRun, romDir)
    # extract the forcing size
    forcingSize = extractForcingSizeFromInputFile(dryRun, romDir)

    if not dryRun:
      # find how many sets of runs are done
      snapsFiles = [i for i in glob.glob("./snaps_vp*")]
      # sort based on the ID
      def func(elem): return int(elem.split('_')[-1])
      snapsFiles = sorted(snapsFiles,key=func)

      linkExe(romDir, exeDir, reconstructSeismoExeName)
      for k,isnap in enumerate(snapsFiles):
        args = "./"+reconstructSeismoExeName +\
          " --podmodes="+vpLsv+' binary' +\
          " --romsize="+str(romSize) +\
          " --romsnaps="+isnap +" binary" +\
          " --fsize="+str(forcingSize) +\
          " --outformat=ascii" +\
          " --gridids=" + gidsString+\
          " --outfileappend=" + str(k)

        print("\nRunning: {} inside {}".format(args, romDir))
        logfile = 'out_reconstruct_seismo_'+str(k)+'.log'
        os.system(args + ' > ' + logfile)
        ## check that it worked
        assert('success' in open(logfile).read().split())

    else:
      print("dryrun:rom: {}".format(romDir))

    os.chdir(owd)



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

  parser.add_argument("-numsamples", "--numsamples",
                      dest="numSamples", type=int,
                      help="Num of samples. Required!")

  parser.add_argument("-mesh", "--mesh",
                      nargs="*", dest="meshSize", type=int,
                      help="Num of points to use along r and theta. Required!")

  #------------------------------
  # parse all args
  #------------------------------
  args = parser.parse_args()
  assert(args.workDir != "empty")
  assert(args.exeDir != "empty")
  workDir  = args.workDir
  drun     = args.dryRun
  numSamples = args.numSamples
  assert( len(list(args.meshSize))==2 )

  # check if working dir exists, if not, make it
  if not os.path.exists(workDir): os.system('mkdir -p ' + workDir)

  # check that dir with exe exits
  if not os.path.exists(args.exeDir):
    sys.exit("The exedir {} provided is empty".format(args.exeDir))

  #------------------------------------------
  # list of energies for computing ROMs
  energies = [99.99999999]
  #energies = [99., 99.999, 99.9995, 99.9999, 99.9999999, 99.99999999]
  #energies = [99.99999, 99.99999999]

  # the parameter's range to use
  pRange     = [31., 69.]

  # this is the size of the state/forcing used when we do the rank2 case
  # becasuse we cannot just set it to however many samples we need since
  # that is not allowed by memory or other constraints.
  # Also, for fom this is small number, but for ROM this can be large.
  maxForcingSizeForRank2Fom = 4
  maxForcingSizeForRank2Rom = 64

  # generate samples of forcing period
  np.random.seed(5231276)
  samples = np.random.uniform(low=float(pRange[0]), high=float(pRange[1]), size=(numSamples))
  print("Sampling values = ", samples)

  meshFullPath = generateMesh(workDir, args.meshSize)
  baseDic = createBaseDicProductionRun(meshFullPath)

  # # **** fom runs, use rank-2 mode ****
  # print("\nDoing fom test rank-2 runs")
  # doFomRunsRank2(drun, workDir, baseDic, args.exeDir, samples,
  #                maxForcingSizeForRank2Fom, trainOrTest="test", runIDStart=0)

  # # # **** training runs for POD ****
  # # # train runs are done at the two bounding values
  # # print("\nDoing rom training runs")
  # # doFomRunsRank1(drun, workDir, baseDic, args.exeDir, pRange,
  # #                trainOrTest="train", runIDStart=0)

  # # **** run or load POD ****
  print("\nDoing svd")
  binaryFlag     = '1' if baseDic['io']['snapshotMatrix']['binary'] else '0'
  podDirFullPath = '/home/fnrizzi/waveRuns/romAccuracyScenario2/train2/pod'
  #podDirFullPath = runSVD(drun, len(pRange), workDir, args.exeDir, binaryFlag)

  # compute list of modes corresponding to target energies
  # for dryRun, just put dummy values 1,2 to test the workflow but nothing is run anyway
  modes = getListOfModesToRun(drun, podDirFullPath, energies) if not drun else [1, 2]
  print('List of modes {}'.format(modes))

  # **** ROM runs ****
  print("\nDoing rom")
  # we can do a single rom run with multiforcing
  doRankTwoRoms(drun, workDir, baseDic, args.exeDir, samples,
                maxForcingSizeForRank2Rom, podDirFullPath, modes)

  # **** reconstruct seismogram from rom data ****
  print("\nReconstructing seismogram from rom")
  reconstructSeismogramRom(drun, workDir, args.exeDir, samples, podDirFullPath)
