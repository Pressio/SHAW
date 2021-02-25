#!/usr/bin/env python

import yaml
import numpy as np
import sys, re, os
import matplotlib.pyplot as plt
from numpy import linalg as la
from matplotlib import cm

def loadSeismogram():
  return np.loadtxt('seismogram_0')

def doPlot(panelId, t, data, angle, depths):
  plt.subplot(panelId)
  plt.grid('on')

  rowMapper = {'5':0, '30':1, '55':2, '80':3}
  row = rowMapper[angle]
  d0 = data['d0'][row, :]
  d1 = data['d1'][row, :]
  d2 = data['d2'][row, :]
  d3 = data['d3'][row, :]

  plt.title("Seismogram for receiver at " + angle+'\u00b0',
            fontsize=15, color='w')

  plt.plot(t, d0, '-o', color='m',
           markerfacecolor='none',
           markersize=1, linewidth=1.8,
           label='With source depth='+depths[0]+' km')

  plt.plot(t, d1, '-s', color='c',
           markerfacecolor='none',
           markersize=1, linewidth=1.8,
           label='With source depth='+depths[1]+' km')

  plt.plot(t, d2, '-*', color='r',
           markerfacecolor='none',
           markersize=1, linewidth=1.8,
           label='With source depth='+depths[2]+' km')

  plt.plot(t, d3, '-v', color='y',
           markerfacecolor='none',
           markersize=1, linewidth=1.8,
           label='With source depth='+depths[3]+' km')

  lg = plt.legend(loc="upper right",
             ncol=1, fontsize=12, labelspacing=.3,
             handletextpad=0.2,
             frameon=False, markerscale=0.75)
  plt.setp(lg.get_texts(), color='w')

  plt.xlim([-50, 2050])
  plt.xticks(np.linspace(0, 2000, 6), color='w')
  plt.ylim([-2.5e-6, 2.5e-6])

  ylab = r'$v_{\phi}(t)$'
  plt.ylabel(ylab, fontsize=15)
  plt.xlabel(r'Time (seconds)', fontsize=15)

  ax = plt.gca()
  mycolor = 'w'
  ax.xaxis.label.set_color(mycolor);
  ax.tick_params(axis='x', colors=mycolor)
  ax.yaxis.label.set_color(mycolor);
  ax.tick_params(axis='y', colors=mycolor)


###############################
if __name__== "__main__":
###############################
  # get yaml file that we use to extract things
  inputs = yaml.safe_load(open("input.yaml"))

  finalTime = float(inputs["general"]["finalTime"])
  dt        = float(inputs["general"]["dt"])
  stations  = [str(i) for i in inputs["io"]["seismogram"]["receivers"]]
  seismoFreq= int(inputs["io"]["seismogram"]["freq"])
  depths    = [str(i) for i in inputs["source"]["signal"]["depth"]]
  fSize     = int(inputs["source"]["signal"]["forcingSize"])

  # set the time axis
  # time = samplingFrequency*dt*numSamples
  # see input.yaml for samplingFrequency and dt
  t = seismoFreq*dt*np.arange(finalTime)


  # get all seismogram data
  # data is an array such that:
  # rows = numOfStations * fSize
  # cols = values at each time sampled
  data0 = loadSeismogram()
  # make a dic for each "batch" of data
  numStations = len(stations)
  data = {'d0' : data0[0:numStations, :],
          'd1' : data0[numStations:2*numStations, :],
          'd2' : data0[2*numStations:3*numStations, :],
          'd3' : data0[3*numStations:, :]}

  f = plt.figure(figsize=(13,10))

  # for each target receiver, plot all realizations
  doPlot(311, t, data, stations[0], depths)
  doPlot(312, t, data, stations[1], depths)
  doPlot(313, t, data, stations[2], depths)
  plt.tight_layout()
  f.savefig('seismogram.png', format="png",
            bbox_inches='tight', dpi=300, transparent=True)

  plt.show()
