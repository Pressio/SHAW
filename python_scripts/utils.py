#!/usr/bin/env python

import sys, os, time, shutil, subprocess
import numpy as np
from argparse import ArgumentParser
from constants import *
from regexes import *

def str2bool(v):
  if isinstance(v, bool):
    return v
  if v.lower() in ('yes', 'true', 't', 'y', '1'):
    return True
  elif v.lower() in ('no', 'false', 'f', 'n', '0'):
    return False
  else:
    raise argparse.ArgumentTypeError('Boolean value expected.')

#=======================================================
def linkExe(workDir, exeDir, exeName):
  # check if exeDir exists and contains the executables
  if not os.path.exists(exeDir):
    print("exeDir {} does not exist".format(exeDir))
    sys.exit(1)

  # check if exe exist
  sourceExe  = exeDir+"/"+exeName
  if not os.path.exists(sourceExe):
    print("exe {} does not exist".format(sourceExe))
    sys.exit(1)

  # delete if already in target dir, and relink
  destExe  = workDir+"/"+exeName
  if os.path.exists(destExe):
    os.system("rm -f {}".format(destExe))

  # relink exes
  os.system("ln -s {} {}".format(sourceExe, destExe))

#=======================================================
def generateMeshIfNeeded(workDir, nr, nth):
  meshDir = workDir + "/mesh" + str(nr) + "x" + str(nth)

  if os.path.exists(meshDir):
    # if mesh exists, do nothing
    print('Mesh {} already exists'.format(meshDir))
  else:
    # if not, generate
    print('Generating mesh {}'.format(meshDir))
    # call script
    owd = os.getcwd()
    os.chdir(meshExeDir)
    args = ("python", meshExe,
            "-nr", str(nr), "-nth", str(nth),
            "-working-dir", workDir)

    popen  = subprocess.Popen(args, stdout=subprocess.PIPE); popen.wait()
    output = popen.stdout.read(); #print("\n", output)
    os.chdir(owd)

#=======================================================
def extractRankFromSvdLogFile(logFilePath, dof):
  # \d{1,} match one or more (any) digits before the .
  # \d{9,} match nine or more (any) digits
  svdRankRegExp = re.compile(r'svd_matrix_'+dof+'_rank = \d{1,}')

  file1 = open(logFilePath, 'r')
  strings = re.search(svdRankRegExp, file1.read())
  file1.close()
  return int(strings.group().split()[2])

#=======================================================
def runExe(workDir, exeDir, exeName, nThreads, forceSpread=False):
  # link exeName from exeDir to workDir
  linkExe(workDir, exeDir, exeName)

  my_env = os.environ.copy()
  if (nThreads ==1):
    my_env["OMP_NUM_THREADS"] = str(1)
    my_env["OMP_PLACES"]="{0}"
    my_env["OMP_PROC_BIND"]="true"

  elif (exeName == fomExeName):
    my_env["OMP_NUM_THREADS"] = str(nThreads)
    my_env["OMP_PLACES"]      ="threads"
    my_env["OMP_PROC_BIND"]   ="spread"

  elif (exeName == romExeName):
    my_env["OMP_NUM_THREADS"] = str(nThreads)
    my_env["OMP_PLACES"]      ="cores"
    my_env["OMP_PROC_BIND"]   ="true"

  if (forceSpread):
    my_env["OMP_PROC_BIND"]="spread"

  args = ("./"+exeName, "input.yaml")

  # launch subprocess
  print("Running: {} inside {}".format(args, workDir))
  logfile = open('out.log', 'w')
  p = subprocess.Popen(args, stdout=logfile, stderr=logfile, env=my_env)
  p.wait()
  logfile.close()
  assert(p.returncode==0)

  # popen = subprocess.Popen(args, stdout=subprocess.PIPE); popen.wait()
  # output = popen.stdout.read(); print("\n", output)
  # # print log file
  # with open("out.log", "w") as text_file:
  #   text_file.write("{}".format(output.splitlines()))

#=======================================================
def is_text(fileIn):
    msg = subprocess.Popen(["file", fileIn], stdout=subprocess.PIPE).communicate()[0]
    return re.search('text', str(msg)) != None
