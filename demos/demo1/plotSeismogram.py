#!/usr/bin/env python

import yaml
import numpy as np
import matplotlib.pyplot as plt

def doPlot(panelId, t, data, key):
  plt.subplot(panelId)
  plt.grid('on')

  plt.plot(t, data[key], '-o', color='r',
             markerfacecolor='none',
             markersize=0, linewidth=2, 
             label="Receiver at " + key+'\u00b0')

  lg = plt.legend(loc="upper right",
             ncol=1, fontsize=15, labelspacing=.3,
             handletextpad=0.2,
             frameon=False, markerscale=0.75)
  plt.setp(lg.get_texts(), color='gray')

  plt.xlim([-50, 2050])
  plt.xticks(np.linspace(0, 2000, 6), color='gray')
  plt.ylim([-1.6e-6, 1.6e-6])

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
  # get yaml file that we use to extract things
  inputs = yaml.safe_load(open("input.yaml"))

  finalTime = float(inputs["general"]["finalTime"])
  dt        = float(inputs["general"]["dt"])
  seismoFreq= int(inputs["io"]["seismogram"]["freq"])
  stations  = [str(i) for i in inputs["io"]["seismogram"]["receivers"]]

  # set the time axis
  # time = samplingFrequency*dt*numSamples
  # see input.yaml for samplingFrequency and dt
  t = seismoFreq*dt*np.arange(finalTime)

  # get data
  data = np.loadtxt('seismogram_0')

  D = {}
  D[stations[0]] = data[0, :]
  D[stations[1]] = data[1, :]
  D[stations[2]] = data[2, :]

  # # set the time axis
  # # time = samplingFrequency*dt*numSamples
  # # see input.yaml for samplingFrequency and dt
  # t = 4*0.25*np.arange(2000)

  f = plt.figure(figsize=(13,10))
  doPlot(311, t, D, stations[0])
  doPlot(312, t, D, stations[1])
  doPlot(313, t, D, stations[2])
  plt.tight_layout()
  f.savefig('seismogram.png', format="png", bbox_inches='tight', dpi=300, transparent=True)

  plt.show()
