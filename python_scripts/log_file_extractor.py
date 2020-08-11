#!/usr/bin/env python

import sys, os, time
import numpy as np
import shutil, subprocess
from constants import *

# impl
def _extractMemBw(logFilePath, key):
  reg = re.compile(r''+key+'Bandwidth(\D+)\d+.\d+')
  file1 = open(logFilePath, 'r')
  strings = re.search(reg, file1.read())
  file1.close()
  assert(strings)
  return np.float64(strings.group().split()[2])

def _extractGflops(logFilePath, key):
  reg = re.compile(r''+key+'GFlop(\D+)\d+.\d+')
  file1 = open(logFilePath, 'r')
  strings = re.search(reg, file1.read())
  file1.close()
  assert(strings)
  return np.float64(strings.group().split()[2])

def _extractPerfTime(logFilePath, key):
  reg = re.compile(r''+key+'Time(\D+)\d+.\d+')
  file1 = open(logFilePath, 'r')
  strings = re.search(reg, file1.read())
  file1.close()
  assert(strings)
  return np.float64(strings.group().split()[2])

# usable
def extractNumDofs(logFilePath, dof):
  reg = re.compile(r'numGpt'+dof+' = \d+')
  file1 = open(logFilePath, 'r')
  strings = re.search(reg, file1.read())
  file1.close()
  assert(strings)
  return int(strings.group().split()[2])

def extractNumTotalDofs(logFilePath):
  nVp = extractNumDofs(logFilePath, 'Vp')
  nSp = extractNumDofs(logFilePath, 'Sp')
  return nVp+nSp

def extractComputIntensity(logFilePath):
  reg = re.compile(r'flops/bytes(\D+)\d+.\d+')
  file1 = open(logFilePath, 'r')
  strings = re.search(reg, file1.read())
  file1.close()
  assert(strings)
  return np.float64(strings.group().split()[2])

def extractMemBw(logFilePath):
  aveV = _extractMemBw(logFilePath, "ave")
  minV = _extractMemBw(logFilePath, "min")
  maxV = _extractMemBw(logFilePath, "max")
  return np.array([aveV, minV, maxV])

def extractGflops(logFilePath):
  aveV = _extractGflops(logFilePath, "ave")
  minV = _extractGflops(logFilePath, "min")
  maxV = _extractGflops(logFilePath, "max")
  return np.array([aveV, minV, maxV])

def extractPerfTimes(logFilePath):
  aveV = _extractPerfTime(logFilePath, "ave")
  minV = _extractPerfTime(logFilePath, "min")
  maxV = _extractPerfTime(logFilePath, "max")
  return np.array([aveV, minV, maxV])

def extractLoopTime(logFilePath):
  reg = re.compile(r'loopTime(\D+)\d+.\d+')
  file1 = open(logFilePath, 'r')
  strings = re.search(reg, file1.read())
  file1.close()
  assert(strings)
  return np.float64(strings.group().split()[2])

def extractDataIoTime(logFilePath):
  reg = re.compile(r'dataCollectionTime(\D+)\d+.\d+')
  file1 = open(logFilePath, 'r')
  strings = re.search(reg, file1.read())
  file1.close()
  assert(strings)
  return np.float64(strings.group().split()[2])

def extractRomJacobianTime(logFilePath):
  reg = re.compile(r'romJac time(\D+)\d+.\d+')
  file1 = open(logFilePath, 'r')
  strings = re.search(reg, file1.read())
  file1.close()
  assert(strings)
  return np.float64(strings.group().split()[2])

def extractFinalTime(logFilePath):
  reg = re.compile(r'finalT(\D+)\d+')
  file1 = open(logFilePath, 'r')
  strings = re.search(reg, file1.read())
  file1.close()
  assert(strings)
  return np.float64(strings.group().split()[2])

def extractNumSteps(logFilePath):
  reg = re.compile(r'numSteps(\D+)\d+')
  file1 = open(logFilePath, 'r')
  strings = re.search(reg, file1.read())
  file1.close()
  assert(strings)
  return int(strings.group().split()[2])

def extractRankTwoFlag(logFilePath):
  reg = re.compile(r'enableRankTwoF_.*')
  file1 = open(logFilePath, 'r')
  strings = re.search(reg, file1.read())
  file1.close()
  assert(strings)
  foundStr = strings.group().split()[2]
  if (foundStr in ["true", "True"]):
    return True
  else:
    return False
