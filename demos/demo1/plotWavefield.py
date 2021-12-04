#!/usr/bin/env python

import numpy as np
import sys, re, os
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
from argparse import ArgumentParser
from numpy import linalg as la

# radius of the earth
earthRadius = 6371. #km
# thickness of the mantle
mantleThickness = 2891. # km
# the radius of the core-mantle boundary
cmbRadius = earthRadius - mantleThickness

def plotCMB(ax):
  # trace the CMB
  cmbTh = np.linspace(0, 2*np.pi, 100)
  cmbRa = cmbRadius*np.ones(100)
  ax.plot(cmbTh, cmbRa, c='k', linestyle='-', linewidth=0.5, zorder=2)

def plotEarthSurf(ax):
  # trace the earth surface
  surfTh = np.linspace(0, 2*np.pi, 100)
  surfRa = earthRadius*np.ones(100)
  ax.plot(surfTh, surfRa, c='k', linewidth=0.5)

#=========================================
def doPlot(th, r, z, figID, bd, outName, title, plotSource=False):
  cm1 = plt.cm.get_cmap('cividis')

  fig1 = plt.figure(figID)
  ax1 = fig1.add_subplot(111, projection='polar')

  h1=ax1.pcolormesh(th, r, z, cmap=cm1, shading = "flat",
                    vmin=bd[0], vmax=bd[1], zorder=1)
  ax1.set_ylim([cmbRadius, earthRadius])
  ax1.set_yticks([]) #[3480, 5701, 6371])
  #plt.yticks(fontsize=13)
  ax1.set_thetamin(-90)
  ax1.set_thetamax(90)
  ax1.set_xticks(np.pi/180. * np.linspace(-90, 90., 7, endpoint=True))
  ax1.set_xticklabels([r'$\pi$', r'$5\pi/6$', r'$4\pi/6$',
                       r'$\pi/2$', r'$2\pi/6$', r'$\pi/6$', r'$0$'],
                      fontsize=11)

  ax1.set_title(title, fontsize=15, color='w')
  ax1.set_rorigin(-1)
  plotEarthSurf(ax1)
  plotCMB(ax1)

  # put a red circle to indicate where to source is located
  # the depth is extracted from the input.yaml
  if plotSource:
    sourceRadius = earthRadius-640. #[km]
    c = ax1.scatter(np.pi/2.01, sourceRadius, c='r', s=15)  
    ax1.text(np.pi/2.01, sourceRadius, "Source", horizontalalignment='center', verticalalignment='top', color='w')  

  mycolor = 'w'
  ax1.xaxis.label.set_color(mycolor);
  ax1.tick_params(axis='x', colors=mycolor)
  ax1.yaxis.label.set_color(mycolor);
  ax1.tick_params(axis='y', colors=mycolor)

  plt.tight_layout()
  fig1.savefig(outName, format="png",bbox_inches='tight', dpi=300, transparent=True)
  plt.show()

###############################
if __name__== "__main__":
###############################
  nr, nth = 200, 1000
  cc = np.loadtxt("./coords_vp.txt")
  th, r = -cc[:,0]+np.pi/2., cc[:, 1]/1000. #m to km
  th, r = th.reshape((nr,nth)), r.reshape((nr,nth))

  fomFile = './state_timestep_1000_vp'
  fomState = np.loadtxt(fomFile, skiprows=1)
  doPlot(th, r, fomState.reshape((nr, nth)), 0,[-5e-9, 5e-9], "wavefield_1000.png", title="t=250 (s)", plotSource=True)
  fomFile = './state_timestep_4000_vp'
  fomState = np.loadtxt(fomFile, skiprows=1)
  doPlot(th, r, fomState.reshape((nr, nth)), 0,[-5e-9, 5e-9], "wavefield_4000.png", title="t=1000 (s)")
  fomFile = './state_timestep_8000_vp'
  fomState = np.loadtxt(fomFile, skiprows=1)
  doPlot(th, r, fomState.reshape((nr, nth)), 0,[-5e-9, 5e-9], "wavefield_8000.png", title="t=2000 (s)")

