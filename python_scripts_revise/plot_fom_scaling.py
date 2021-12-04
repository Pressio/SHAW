#!/usr/bin/env python

import sys, os, time, yaml
import pprint as pp
import subprocess, math
import numpy as np
import os.path, re
from scipy import stats
from argparse import ArgumentParser
import matplotlib.pyplot as plt

np.set_printoptions(edgeitems=10, linewidth=100000)

# peak BW for mrstem (GB) = 21333.33 * 4 at time of doing these runs
peakMemBW = 85.

myAlpha = .9

# num of dofs for each case for plotting
#785152, 3143168, 12577792, 50321408
meshLabelsPlot = [r'$0.78~\cdot~10^6$', 
                  r'$3~\cdot~10^6$', 
                  r'$12~\cdot~10^6$', 
                  r'$50~\cdot~10^6$']

nThreads = [2, 8, 36]
colors = {2:'#009286', 8:'#ff9e11', 36:'#cd4d84'}

fSizes = [1,4,16,48]

#=====================================================================
def computeMetricValue(lineData, currValueF, metric, stat):
  if metric == "mem":
    if stat == "ave": return lineData[4]
    elif stat=="min": return lineData[5]
    elif stat=="max": return lineData[6]

  elif metric == "cpu":
    if stat == "ave": return lineData[7]
    elif stat=="min": return lineData[8]
    elif stat=="max": return lineData[9]

  elif metric == "itertime":
    if stat == "ave": return lineData[10]
    elif stat=="min": return lineData[11]
    elif stat=="max": return lineData[12]

  elif metric == "looptime":
    return lineData[13]

#=====================================================================
def createDataDic(data, metric, stat):
  all = {}
  for nt in nThreads:
    dic = {}
    for i in range(data.shape[0]):
      # number of threads and number of modes
      thisNumThr   = int(data[i][0])
      thisValF     = int(data[i][1])

      if thisNumThr == nt and thisValF in fSizes:
          value = computeMetricValue(data[i,:], thisValF, metric, stat)
          if thisValF in dic: dic[thisValF].append(value)
          else: dic[thisValF] = [value]
    all[nt] = dic

  return all

#=====================================================================
def plotBarSet(ax, xLoc, width, f, dic, myColor):
  val = dic[f]
  ax.bar(xLoc, val, width, alpha=myAlpha, color=myColor, edgecolor='none', zorder=5)

#=====================================================================
def plotBar(dataDic, meshLabels, nThreads, metric, stat):
  # number of mesh sizes to deal with
  numMeshes = len(meshLabels)
  # Setting the positions and width for the bars
  posArray = range(numMeshes)
  pos = list(posArray)

  width = 0.45
  plt.rc('axes', axisbelow=True)

  fig, ax = plt.subplots(figsize=(9,6))
  plt.grid()
  ax2 = ax.twiny()
  fig.subplots_adjust(bottom=0.25)

  gigi = [0.25, 6.5, 12.75, 19.]

  xTicksBars, xTlabels = [], []
  count=0
  for k,v in dataDic.items():
    for i,f in enumerate(fSizes):
      #x locations for the bars
      shift = width*i*3.5

      xLoc = [p+shift+0.455*count+gigi[k] for k,p in enumerate(pos)]

      plotBarSet(ax, xLoc, width, f, v, colors[k])
      xTicksBars += [p+shift+0.475+gigi[k] for k,p in enumerate(pos)]
      xTlabels += [str(f) for i in range(numMeshes)]
    count+=1

  for nt in nThreads:
    ax.bar(100, 1, width, alpha=myAlpha, color=colors[nt],
           edgecolor='none', zorder=-1, label='threads='+str(nt))

  l = ax.legend(loc="upper center", ncol=5, fontsize=13, frameon=False)
  for text in l.get_texts():
    text.set_color("gray")

  # remove the vertical lines of the grid
  ax.xaxis.grid(which="major", color='None', linestyle='-.', linewidth=0, zorder=0)

  ax.xaxis.set_ticks_position('bottom')
  ax.xaxis.set_label_position('bottom')
  ax.set_xticks(xTicksBars)
  ax.set_xticklabels(xTlabels, fontsize=15, color='gray')
  ax.xaxis.set_tick_params(rotation=0)

  ax.set_xlabel('Number of simultaneous trajectories (M)', fontsize=16, color='gray')
  ax.set_xlim(min(pos)-0.2, max(pos)+width*56)

  if metric =="mem":
    ax.set_yscale('log')
    ax.set_ylabel("Memory Bandwith (GB/s)", fontsize=18)
    ax.set_ylim([1e-1, 1000])
    ax.set_yticks([1e-1, 1, 10, 100, 1000])
    ax.tick_params(axis='y', which='major', labelsize=15, color='gray')
    ax.tick_params(axis='y', which='minor', labelsize=13, color='gray')

    # # plot peak theoretical mem BW
    # ax.plot([min(pos)-0.2, max(pos)+width*70],
    #         [peakMemBW, peakMemBW], '--k', linewidth=1.2, zorder=7)
    # ax.text((min(pos)+width+max(pos)+width*75)*0.45,
    #         peakMemBW+12, 'Machine\'s theoretical peak', fontsize=15)

  elif metric=='cpu':
    ax.set_yscale('log')
    ax.set_ylabel("GFlops", fontsize=18, color='gray')
    ax.set_ylim([1e-1, 1e4])
    ax.set_yticks([1e-1, 1, 10, 1e2, 1e3, 1e4])
    ax.tick_params(axis='y', which='major', labelsize=15, color='gray')
    ax.tick_params(axis='y', which='minor', labelsize=13, color='gray')

  elif metric =="itertime":
    ax.set_yscale('log')
    ax.set_ylim([1e-1, 1e4])
    ax.tick_params(axis='y', which='major', labelsize=15, color='gray')
    ax.tick_params(axis='y', which='minor', labelsize=13, color='gray')
    if stat == 'ave': pref = 'Average'
    elif stat=='min': pref = 'Min'
    elif stat=='max': pref = 'Max'
    ax.set_ylabel(pref+" time (ms)/timestep", fontsize=18, color='gray')


  # ticks for the meshes
  meshTicks = [3.5, 11., 18.75, 26.35]
  ax2.set_xticks(meshTicks)
  ax2.xaxis.set_ticks_position('bottom')
  ax2.xaxis.set_label_position('bottom')
  ax2.spines['bottom'].set_position(('outward', 65))
  ax2.set_xlabel('Total degrees of freedom (N)', fontsize=16, color='gray')
  ax2.set_xticklabels(meshLabels, fontsize=16, color='gray')
  ax2.set_xlim(min(pos), max(pos)+width*60)
  ax2.set_axisbelow(True)

  ax.tick_params(axis='y', colors='gray')

  ax.xaxis.label.set_color('gray')
  ax.yaxis.label.set_color('gray')
  ax2.xaxis.label.set_color('gray')
  ax2.yaxis.label.set_color('gray')

  plt.tight_layout()
  fileName = "fom_"+metric+"_"+stat+".png"
  fig.savefig('./plots/'+fileName, format="png", bbox_inches='tight', 
    dpi=300, transparent=True)
  plt.show()

#=====================================================================
def main(dataFile, metric, stat):
  data = np.loadtxt(dataFile)
  dataDic = createDataDic(data, metric, stat)
  #print(dataDic)
  pp.pprint(dataDic)

  plotBar(dataDic, meshLabelsPlot, nThreads, metric, stat)
  plt.show()


#////////////////////////////////////////////
if __name__== "__main__":
#////////////////////////////////////////////
  parser = ArgumentParser()
  # parser.add_argument("-file", "--file",
  #                     dest="dataFile",
  #                     help="where to get data from\n")

  parser.add_argument("-metric", "--metric",
                      dest="metric", default="mem",
                      help="Choices: mem, cpu, itertime \n")

  parser.add_argument("-stat", "--stat",
                      dest="stat", default="ave",
                      help="ave, min or max\n")

  args = parser.parse_args()

  assert(args.metric in ['mem', 'cpu', 'itertime'])
  main('./data/fom_scaling_final.txt', args.metric, args.stat)

#////////////////////////////////////////////
