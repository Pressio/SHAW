#!/usr/bin/env python

import matplotlib.pyplot as plt
import sys, os, time
import numpy as np
from numpy import linspace, meshgrid
from matplotlib import cm

from earth_params import earthRadius, mantleThickness, cmbRadius

#----------------------------------------------------------------------------------------
def plotCMB(ax):
  # trace the CMB
  cmbTh = np.linspace(0, 2*np.pi, 100)
  cmbRa = cmbRadius*np.ones(100)
  ax.plot(cmbTh, cmbRa, c='b', linestyle='--', linewidth=0.5)

def plotEarthSurf(ax):
  # trace the earth surface
  surfTh = np.linspace(0, 2*np.pi, 100)
  surfRa = earthRadius*np.ones(100)
  ax.plot(surfTh, surfRa, c='b', linewidth=0.5)
