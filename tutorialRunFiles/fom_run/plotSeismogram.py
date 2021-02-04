#!/usr/bin/env python

import numpy as np
import sys, re, os
import matplotlib.pyplot as plt
from numpy import linalg as la
from matplotlib import cm

# this script is supposed to be used
# as described in the tutorial ../docs/run.md


def loadSeismogram():
  return np.loadtxt('seismogram_0')

def doPlot(panelId, t, data, key):
  plt.subplot(panelId)
  plt.grid('on')

  plt.plot(t, data[key], '-o', color='k',
             markerfacecolor='none',
             markersize=3, linewidth=1, label=key+'\u00b0')

  plt.legend(loc="upper right",
             ncol=1, fontsize=15, labelspacing=.3,
             handletextpad=0.2,
             frameon=False, markerscale=0.75)

  plt.xlim([-50, 2050])
  plt.xticks(np.linspace(0, 2000, 6))
  plt.ylim([-1.6e-6, 1.6e-6])

  ylab = r'$v_{\phi}(t)$'
  plt.ylabel(ylab, fontsize=15)
  plt.xlabel(r'Time (seconds)', fontsize=15)


###############################
if __name__== "__main__":
###############################
  # get data
  data = loadSeismogram()

  # get data for each receiver:
  # inside the input.yaml file, we set receivers
  # at 5,30,55,80,105, etc degrees, but here only show a few
  D = {}
  D['5']  = data[0, :]
  D['30'] = data[1, :]
  D['55'] = data[2, :]

  # set the time axis
  # time = samplingFrequency*dt*numSamples
  # see input.yaml for samplingFrequency and dt
  t = 4*0.25*np.arange(2000)

  f = plt.figure(figsize=(13,10))
  doPlot(311, t, D, '5')
  doPlot(312, t, D, '30')
  doPlot(313, t, D, '80')
  plt.tight_layout()
  f.savefig('seismogram.png', format="png", bbox_inches='tight', dpi=300)

  plt.show()
