#!/usr/bin/env python

from random import randrange, uniform
import re, sys, os, time, yaml, copy
import numpy as np
from argparse import ArgumentParser
import matplotlib.pyplot as plt

np.set_printoptions(edgeitems=10, linewidth=100000)

sizeofdouble = 8 # bytes

#--------------------------------------------------------------
def memoryFootprintFom(vpStateSize, spStateSize, forcingSize,
                       numSteps, samplingFreq):
  # neglect Jacobians for now

  statesBytes  = (vpStateSize + spStateSize)*sizeofdouble

  # forcing we only store the time series
  forcingBytes = sizeofdouble * numSteps * forcingSize

  # observer
  numObs = numSteps/samplingFreq
  obsVpBytes = sizeofdouble * vpStateSize * numObs * forcingSize
  obsSpBytes = sizeofdouble * spStateSize * numObs * forcingSize

  return statesBytes+forcingBytes+obsVpBytes+obsSpBytes


#--------------------------------------------------------------
def memoryFootprintRom(romSize, vpStateSize, spStateSize,
                       forcingSize, numSteps, samplingFreq):
  # Jacobians
  jacsBytes = 2.*romSize*romSize*sizeofdouble

  statesBytes  = (2.*romSize)*sizeofdouble

  # forcing we store the time series and also the forcing matrix
  forcingBytes = sizeofdouble*romSize*forcingSize + sizeofdouble * numSteps * forcingSize

  # observer
  numObs = numSteps/samplingFreq
  obsBytes = 2. * sizeofdouble * romSize * numObs * forcingSize

  return jacsBytes+obsBytes+statesBytes+forcingBytes



dt = 0.1 # seconds
simulationTime = 2000. # seconds
samplingFreq = 20 #num of steps
totalSteps = int(simulationTime/dt)
forcingSize = 8.

fomVpSize = 1048576
fomSpSize = 2094592

romSize = 436

M1 = memoryFootprintFom(fomVpSize, fomSpSize, forcingSize, totalSteps, samplingFreq)
print("fom (GB) = {}".format(M1/(1024.*1024.*1024.)))

M2 = memoryFootprintRom(romSize, fomVpSize, fomSpSize, forcingSize, totalSteps, samplingFreq)
print("rom (GB) = {}".format(M2/(1024.*1024.*1024.)))
