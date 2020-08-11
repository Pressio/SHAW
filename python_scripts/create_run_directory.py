#!/usr/bin/env python

import sys, os, time
import numpy as np
import shutil, re

def getMeshString(D):
  return D['general']['meshDir'].split('/')[-1]

def getDtString(D):
  return str(D['general']['dt'])

def getTotalTimeString(D):
  return str(D['general']['finalTime'])

def getSnapsOn(D):
  hasIo = True if 'io' in D.keys() else False
  hasSM = False
  if hasIo:
    if 'snapshotMatrix' in D['io'].keys(): hasSM = True
  return 'true' if hasIo and hasSM else 'false'

def getSeismoOn(D):
  hasIo = True if 'io' in D.keys() else False
  hasSE = False
  if hasIo:
    if 'seismogram' in D['io'].keys(): hasSE = True
  return 'true' if hasIo and hasSE else 'false'

def getForcingRank(D):
  if 'sampling' not in D.keys():
    return str(1)
  else:
    fSize = D['sampling']['forcingSize']
    return str(fSize)

def getMaterialType(D):
  return D['material']['kind']

def getNumModes(D, dof):
  return D['rom'][dof]['numModes']


def _createRunDirectoryImpl(isFom, nThreads, runDic):
  dd = "fom" if isFom == True else "rom"
  dd += "_"+getMeshString(runDic)
  dd += "_nThreads_"+str(nThreads)
  dd += "_dt_"+getDtString(runDic)
  dd += "_T_"+getTotalTimeString(runDic)
  dd += "_snaps_"+getSnapsOn(runDic)

  # the seismogram is only computed for FOM
  if (isFom):
    dd += "_seismo_"+getSeismoOn(runDic)

  dd += "_mat_"+getMaterialType(runDic)
  dd += "_fRank_"+getForcingRank(runDic)
  return dd

def createFomRunDirectory(workDir, nThreads, runDic, appendable=""):
  subDir = _createRunDirectoryImpl(True, nThreads, runDic)
  if (appendable == ""):
    return workDir+"/"+subDir
  else:
    return workDir+"/"+subDir+appendable

def createRomRunDirectory(workDir, nThreads, D, appendable=""):
  subDir = _createRunDirectoryImpl(False, nThreads, D)
  nModesVp = getNumModes(D, 'velocity')
  nModesSp = getNumModes(D, 'stress')
  subDir += "_nPod_"+str(nModesVp)+"_"+str(nModesSp)
  if (appendable == ""):
    return workDir+"/"+subDir
  else:
    return workDir+"/"+subDir+appendable
