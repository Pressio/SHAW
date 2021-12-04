#!/usr/bin/env python

import numpy as np
import sys, re, os
import matplotlib.pyplot as plt
from numpy import linalg as la
from matplotlib import cm

def loadSeismograms():
  return {'d0' : np.loadtxt('seismogram_0'),
          'd1' : np.loadtxt('seismogram_1'),
          'd2' : np.loadtxt('seismogram_2'),
          'd3' : np.loadtxt('seismogram_3')}

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
            fontsize=15, color='gray')

  plt.plot(t, d0, '-o', color='m',
           markerfacecolor='none',
           markersize=0, linewidth=1.8,
           label='With source depth='+depths[0]+' km')

  plt.plot(t, d1, '-s', color='c',
           markerfacecolor='none',
           markersize=0, linewidth=1.8,
           label='With source depth='+depths[1]+' km')

  plt.plot(t, d2, '-*', color='r',
           markerfacecolor='none',
           markersize=0, linewidth=1.8,
           label='With source depth='+depths[2]+' km')

  plt.plot(t, d3, '-v', color='y',
           markerfacecolor='none',
           markersize=0, linewidth=1.8,
           label='With source depth='+depths[3]+' km')

  lg = plt.legend(loc="upper right",
             ncol=1, fontsize=12, labelspacing=.3,
             handletextpad=0.2,
             frameon=False, markerscale=0.75)
  plt.setp(lg.get_texts(), color='gray')

  plt.xlim([-50, 2050])
  plt.xticks(np.linspace(0, 2000, 6), color='gray')
  plt.ylim([-2.5e-6, 2.5e-6])

  ylab = r'$v_{\phi}(t)$'
  plt.ylabel(ylab, fontsize=15)
  plt.xlabel(r'Time (seconds)', fontsize=15)

  ax = plt.gca()
  mycolor = 'gray'
  ax.xaxis.label.set_color(mycolor);
  ax.tick_params(axis='x', colors=mycolor)
  ax.yaxis.label.set_color(mycolor);
  ax.tick_params(axis='y', colors=mycolor)


###############################
if __name__== "__main__":
###############################
  # get data for all realizations
  data = loadSeismograms()

  depths = ['240', '440', '540', '700']

  # set the time axis
  # time = samplingFrequency*dt*numSamples
  # see input.yaml for samplingFrequency and dt
  t = 4*0.25*np.arange(2000)

  f = plt.figure(figsize=(13,10))

  # for each target receiver, plot all realizations
  doPlot(311, t, data, '5',  depths)
  doPlot(312, t, data, '30', depths)
  doPlot(313, t, data, '55', depths)
  plt.tight_layout()
  f.savefig('seismogram.png', format="png",
            bbox_inches='tight', dpi=300, transparent=True)

  plt.show()
