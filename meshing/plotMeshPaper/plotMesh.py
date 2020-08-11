#!/usr/bin/env python

import matplotlib.pyplot as plt
import sys, os, time
import numpy as np
from numpy import linspace, meshgrid
from matplotlib import cm
import collections
import random
from argparse import ArgumentParser

# radius of the earth
earthRadius = 6371. #km
# thickness of the mantle
mantleThickness = 2891. # km
# the radius of the core-mantle boundary
cmbRadius = earthRadius - mantleThickness

def plotThetaLine(ax):
  ax.plot(np.linspace(np.pi/2., np.pi/4, 100),
          1800*np.ones(100), c='k', linestyle='-', linewidth=1.25)

def plotCMB(ax):
  # trace the CMB
  cmbTh = np.linspace(0, 2*np.pi, 100)
  cmbRa = cmbRadius*np.ones(100)
  ax.plot(cmbTh, cmbRa, c='k', linestyle='-', linewidth=1.25)

def plotEarthSurf(ax):
  # trace the earth surface
  surfTh = np.linspace(0, 2*np.pi, 100)
  surfRa = earthRadius*np.ones(100)
  ax.plot(surfTh, surfRa, c='k', linewidth=1.25)

def doPlot(workDir):
  vpD = np.loadtxt(workDir+'/graph_vp.dat')
  spD = np.loadtxt(workDir+'/graph_sp.dat')
  # convert rad from m to km
  #vpD[:,3]/= 1000.
  #spD[:,4]/= 1000.

  srpD = spD[spD[:,1]==1]
  stpD = spD[spD[:,1]==2]

  fig = plt.figure(0)
  ax = fig.add_subplot(111, frameon=False, projection='polar')
  #ax.set_rlabel_position(90)
  #ax.set_aspect(1.0)

  plotEarthSurf(ax)
  plotCMB(ax)
  plotThetaLine(ax)

  # the minus sign for thetat below is needed to plot things so that
  # it shows the grid on the right half of the circle
  # TODO: fix labels of theta

  v_th, v_r = -vpD[:,2]+np.pi*0.5, vpD[:,3]
  ax.set_ylim([0, earthRadius+200])
  ax.set_yticks(v_r)

  ax.set_xticks(v_th) #np.pi/180. * np.linspace(-90, 120, 7, endpoint=False))
  ax.set_xticklabels([r'$0$','',r'$\pi/4$','',r'$\pi/2$','',r'$3\pi/4$','',r'$\pi$'], fontsize=15)
  #ax.set_yticklabels([])

  #ax.set_rlabel_position(0)
  ax.set_thetamin(-90)
  ax.set_thetamax(90)

  # arr1 = plt.arrow(0, 0, 0, 1000, alpha = 1,
  #                  edgecolor = 'black', facecolor = 'black', lw = 2, zorder = 5, clip_on=False,
  #                  head_width=0.5, head_length=0.5, capstyle='round')
  #ax.quiver((0, 0), (1, 0), color='black', clip_on=False, linewidth=0.01)
  ax.quiver((1, 0), (1,0), color='black', clip_on=False, linewidth=0.00001, minlength=1, scale=7.5)

  plt.text(0.3, 0.5, r'$r$', transform=ax.transAxes, fontsize=20)
  plt.text(0.28, 0.65, r'$\theta$', transform=ax.transAxes, fontsize=20)

  ms=25
  mm=['o','s','<']

  v_th, v_r = -vpD[:,2]+np.pi*0.5, vpD[:,3]
  ax.scatter(v_th, v_r, marker=mm[0], s=25, edgecolor='k', facecolor='red', clip_on=False, zorder=5)

  sr_th, sr_r = -srpD[:,3]+np.pi*0.5, srpD[:,4]
  ax.scatter(sr_th, sr_r, marker=mm[1], s=ms, edgecolor='k', facecolor='gray', clip_on=False, zorder=5)

  st_th, st_r = -stpD[:,3]+np.pi*0.5, stpD[:,4]
  ax.scatter(st_th, st_r, marker=mm[2], s=ms, edgecolor='k', facecolor='gray', zorder=5)

  ax.set_axisbelow(True)
  plt.savefig("mesh.png", bbox_inches='tight', dpi=600)

#     th_fm = -np.asarray([v[0] for (k,v) in coordsSp.items()])+np.pi*0.5
#     r_fm  = np.asarray([v[1] for (k,v) in coordsSp.items()])
#     gids  = np.asarray([k for (k,v) in coordsSp.items()])
#     ax.scatter(th_fm, r_fm, marker='s', s=4, edgecolor='k', facecolor='k')
#     for i in range(len(gids)):
#       ax.text(th_fm[i], r_fm[i], str( gids[i] ),
#               verticalalignment='top', horizontalalignment='center',
#               fontsize=10, color='k')
  plt.show()


###############################
if __name__== "__main__":
###############################
  parser = ArgumentParser()
  parser.add_argument("-wdir", "--wdir", dest="workDir", default=".",
                      help="Where mesh data lives")

  args = parser.parse_args()
  doPlot(args.workDir)
